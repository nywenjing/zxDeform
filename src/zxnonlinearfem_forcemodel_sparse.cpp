#include "zxnonlinearfem_forcemodel_sparse.h"
#include "zxmath.h"
zxNonlinearFEM_ForceModel_Sparse::zxNonlinearFEM_ForceModel_Sparse(zxSolidMesh::Ptr mesh,zxMaterial::Ptr material,real low_clamp,real upper_clamp)
{
    m_mesh = mesh;
    m_material = material;
    m_dim = 3 * m_mesh->get_num_nodes();
    m_low_svd_diag = low_clamp;
    m_upper_svd_diag = upper_clamp;

    std::vector<Eigen::Triplet<real>> k_trips;
    int nnz = 0;
    for(size_t el = 0 ; el < m_mesh->get_num_elements(); el++)
    {
        size_t nn = m_mesh->get_element(el)->get_num_nodes();
        nnz += 3 * nn * 3 * nn;
    }
    k_trips.reserve(nnz);

    for(size_t el = 0 ; el < m_mesh->get_num_elements(); el++)
    {
        zxElement::Ptr element = m_mesh->get_element(el);
        for(int i = 0; i < element->get_num_nodes(); i++)
            for(int j = 0; j < element->get_num_nodes(); j++)
                for(int ii = 0; ii < 3; ii++)
                    for(int jj = 0; jj < 3; jj++)
                    {
                        int r_id = 3 * element->get_node(i)->m_id + ii;
                        int c_id = 3 * element->get_node(j)->m_id + jj;

                        k_trips.push_back(Eigen::Triplet<real>(r_id,c_id,0.0));
                    }
    }
    m_spK_topology.resize(m_dim,m_dim);

    m_spK_topology.setFromTriplets(k_trips.begin(),k_trips.end());
    m_isSparse = true;
}


void zxNonlinearFEM_ForceModel_Sparse::updatePosition(const Eigen::VectorXd &u)
{
    for(size_t i = 0; i < m_mesh->get_num_nodes(); i++)
        for(int j = 0; j < 3; j++)
            m_mesh->get_node(i)->rt[j] = m_mesh->get_node(i)->r0[j] + u[3 * i + j];

    for(size_t el = 0; el < m_mesh->get_num_elements(); el++)
    {
        zxElement::Ptr element = m_mesh->get_element(el);
        for(int g_id = 0; g_id < element->get_num_gaussian_points(); g_id++)
        {
            zxMaterialPoint::Ptr p_mat = element->get_material_point(g_id);

            real* Grn = &element->get_first_derive_r(0,g_id);
            real* Gsn = &element->get_first_derive_s(0,g_id);
            real* Gtn = &element->get_first_derive_t(0,g_id);

            mat3d& Ji = p_mat->m_invJac0;
            mat3d& F = p_mat->m_defgrad;
            F.setZero();
            for(size_t n_id = 0; n_id < element->get_num_nodes(); n_id++)
            {
                double Gri = Grn[n_id];
                double Gsi = Gsn[n_id];
                double Gti = Gtn[n_id];

                double x = element->get_node(n_id)->rt[0];
                double y = element->get_node(n_id)->rt[1];
                double z = element->get_node(n_id)->rt[2];

                double GX = Ji(0,0)*Gri+Ji(1,0)*Gsi+Ji(2,0)*Gti;
                double GY = Ji(0,1)*Gri+Ji(1,1)*Gsi+Ji(2,1)*Gti;
                double GZ = Ji(0,2)*Gri+Ji(1,2)*Gsi+Ji(2,2)*Gti;

                F(0,0) += GX*x; F(0,1) += GY*x; F(0,2) += GZ*x;
                F(1,0) += GX*y; F(1,1) += GY*y; F(1,2) += GZ*y;
                F(2,0) += GX*z; F(2,1) += GY*z; F(2,2) += GZ*z;
            }

            zxModified_SVD(p_mat->m_defgrad,p_mat->m_svd_U,p_mat->m_svd_diag,p_mat->m_svd_V);

            for(size_t ii = 0; ii < 3; ii++)
            {
                p_mat->m_svd_diag[ii] = std::min(p_mat->m_svd_diag[ii],m_upper_svd_diag);
                p_mat->m_svd_diag[ii] = std::max(p_mat->m_svd_diag[ii],m_low_svd_diag);
            }
            zx_vega_ComputeDiagonalPstress(m_material,p_mat->m_svd_diag,p_mat->m_diag_P);

            mat3d P;
            P.setZero();
            for(int ii = 0; ii < 3; ii++)
                P(ii,ii) = p_mat->m_diag_P[ii];
            p_mat->m_Pstress = p_mat->m_svd_U * P * p_mat->m_svd_V.transpose();
        }
    }

}

std::vector<size_t>    zxNonlinearFEM_ForceModel_Sparse::get_element_node_id(int el)
{
    std::vector<size_t> elId(m_mesh->get_element(el)->get_num_nodes());
    for(size_t i = 0; i < m_mesh->get_element(el)->get_num_nodes(); i++)
        elId[i] = m_mesh->get_element(el)->get_node(i)->m_id;

    return elId;

}

void zxNonlinearFEM_ForceModel_Sparse::computeElementForce(int el, Eigen::VectorXd &force, Eigen::VectorXd &u)
{
    zxElement::Ptr element = m_mesh->get_element(el);

    for(size_t i = 0; i < element->get_num_nodes(); i++)
        for(int j = 0; j < 3; j++)
            element->get_node(i)->rt[j] = element->get_node(i)->r0[j] + u[3 * i + j];

    force.resize(element->get_num_nodes() * 3);
    force.setZero();

    for(int g_id = 0; g_id < element->get_num_gaussian_points(); g_id++)
    {
        zxMaterialPoint::Ptr p_mat = element->get_material_point(g_id);

        real* Grn = &element->get_first_derive_r(0,g_id);
        real* Gsn = &element->get_first_derive_s(0,g_id);
        real* Gtn = &element->get_first_derive_t(0,g_id);

        mat3d& Ji = p_mat->m_invJac0;
        mat3d& F = p_mat->m_defgrad;
        F.setZero();
        for(size_t n_id = 0; n_id < element->get_num_nodes(); n_id++)
        {
            double Gri = Grn[n_id];
            double Gsi = Gsn[n_id];
            double Gti = Gtn[n_id];

            double x = element->get_node(n_id)->rt[0];
            double y = element->get_node(n_id)->rt[1];
            double z = element->get_node(n_id)->rt[2];

            double GX = Ji(0,0)*Gri+Ji(1,0)*Gsi+Ji(2,0)*Gti;
            double GY = Ji(0,1)*Gri+Ji(1,1)*Gsi+Ji(2,1)*Gti;
            double GZ = Ji(0,2)*Gri+Ji(1,2)*Gsi+Ji(2,2)*Gti;

            F(0,0) += GX*x; F(0,1) += GY*x; F(0,2) += GZ*x;
            F(1,0) += GX*y; F(1,1) += GY*y; F(1,2) += GZ*y;
            F(2,0) += GX*z; F(2,1) += GY*z; F(2,2) += GZ*z;
        }

        zxModified_SVD(p_mat->m_defgrad,p_mat->m_svd_U,p_mat->m_svd_diag,p_mat->m_svd_V);

        for(size_t ii = 0; ii < 3; ii++)
        {
            p_mat->m_svd_diag[ii] = std::min(p_mat->m_svd_diag[ii],m_upper_svd_diag);
            p_mat->m_svd_diag[ii] = std::max(p_mat->m_svd_diag[ii],m_low_svd_diag);
        }
        zx_vega_ComputeDiagonalPstress(m_material,p_mat->m_svd_diag,p_mat->m_diag_P);

        mat3d P;
        P.setZero();
        for(int ii = 0; ii < 3; ii++)
            P(ii,ii) = p_mat->m_diag_P[ii];
        p_mat->m_Pstress = p_mat->m_svd_U * P * p_mat->m_svd_V.transpose();
    }

    for (size_t g_id =0; g_id <element->get_num_gaussian_points(); ++g_id)
    {
        double Gx, Gy, Gz;

        zxMaterialPoint::Ptr m_pt = element->get_material_point(g_id);

        mat3d& Ji0 = m_pt->m_invJac0;
        double detJ0 = m_pt->m_detJac0;
        mat3d& P = m_pt->m_Pstress;

        detJ0 *= element->get_gaussian_weight(g_id);

        real* Gr = &element->get_first_derive_r(0,g_id);
        real* Gs = &element->get_first_derive_s(0,g_id);
        real* Gt = &element->get_first_derive_t(0,g_id);

        for(size_t n_id = 0; n_id < element->get_num_nodes(); n_id++)
        {

            Gx = Ji0(0,0)*Gr[n_id]+Ji0(1,0)*Gs[n_id]+Ji0(2,0)*Gt[n_id];
            Gy = Ji0(0,1)*Gr[n_id]+Ji0(1,1)*Gs[n_id]+Ji0(2,1)*Gt[n_id];
            Gz = Ji0(0,2)*Gr[n_id]+Ji0(1,2)*Gs[n_id]+Ji0(2,2)*Gt[n_id];

            //size_t gn_id = element->get_node(n_id)->m_id;

            force[3 * n_id + 0] += (  Gx * P(0,0) +
                                       Gy * P(0,1) +
                                       Gz * P(0,2)) * detJ0;
            force[3 * n_id + 1] += (  Gx * P(1,0) +
                                       Gy * P(1,1) +
                                       Gz * P(1,2)) * detJ0;

            force[3 * n_id + 2] += (  Gx * P(2,0) +
                                       Gy * P(2,1) +
                                       Gz * P(2,2)) * detJ0;
        }
    }


}

void zxNonlinearFEM_ForceModel_Sparse::computeForce(Eigen::VectorXd &force)
{
    force.resize(getForceDimension());
    force.setZero();

    //    for(size_t el = 0; el < m_mesh->get_num_elements(); el++)
    //    {
    //        zxElement::Ptr element = m_mesh->get_element(el);

    //        for(size_t g_id = 0; g_id < element->get_num_gaussian_points(); g_id++)
    //        {
    //            zxMaterialPoint::Ptr p_mat = element->get_material_point(g_id);

    //            for(size_t n_id = 0; n_id < element->get_num_nodes(); n_id++)
    //            {
    //                const vec3d r = element->get_node(n_id)->rp;

    //                real dg_drst[3] = {element->get_first_derive_r(n_id,g_id),
    //                                   element->get_first_derive_s(n_id,g_id),
    //                                   element->get_first_derive_t(n_id,g_id)};

    //                for(size_t i = 0; i < 3; i++)
    //                {
    //                    size_t dof_id = 3 * element->get_node(n_id)->m_id + i;

    //                    for(size_t ll = 0; ll < 3; ll++)
    //                        for(size_t kk = 0; kk < 3; kk++)
    //                            force[dof_id] -= p_mat->m_invJac0(kk,ll) * dg_drst[ll] * p_mat->m_Pstress(i,kk) * p_mat->m_detJac0;
    //                }

    //            }
    //        }

    //    }



    double Gx, Gy, Gz;

    for(size_t el = 0; el < m_mesh->get_num_elements(); el++)
    {
        zxElement::Ptr element = m_mesh->get_element(el);

        for (size_t g_id =0; g_id <element->get_num_gaussian_points(); ++g_id)
        {
            zxMaterialPoint::Ptr m_pt = element->get_material_point(g_id);

            mat3d& Ji0 = m_pt->m_invJac0;
            double detJ0 = m_pt->m_detJac0;
            mat3d& P = m_pt->m_Pstress;

            detJ0 *= element->get_gaussian_weight(g_id);

            real* Gr = &element->get_first_derive_r(0,g_id);
            real* Gs = &element->get_first_derive_s(0,g_id);
            real* Gt = &element->get_first_derive_t(0,g_id);

            for(size_t n_id = 0; n_id < element->get_num_nodes(); n_id++)
            {

                Gx = Ji0(0,0)*Gr[n_id]+Ji0(1,0)*Gs[n_id]+Ji0(2,0)*Gt[n_id];
                Gy = Ji0(0,1)*Gr[n_id]+Ji0(1,1)*Gs[n_id]+Ji0(2,1)*Gt[n_id];
                Gz = Ji0(0,2)*Gr[n_id]+Ji0(1,2)*Gs[n_id]+Ji0(2,2)*Gt[n_id];

                size_t gn_id = element->get_node(n_id)->m_id;

                force[3 * gn_id + 0] += (  Gx * P(0,0) +
                                           Gy * P(0,1) +
                                           Gz * P(0,2)) * detJ0;
                force[3 * gn_id + 1] += (  Gx * P(1,0) +
                                           Gy * P(1,1) +
                                           Gz * P(1,2)) * detJ0;

                force[3 * gn_id + 2] += (  Gx * P(2,0) +
                                           Gy * P(2,1) +
                                           Gz * P(2,2)) * detJ0;
            }
        }
    }
}

void zxNonlinearFEM_ForceModel_Sparse::computeTangent(zxSparseMatrix &tangentK)
{
    //    for(size_t el = 0; el < m_mesh->get_num_elements(); el++)
    //    {
    //        zxElement::Ptr element = m_mesh->get_element(el);

    //        Eigen::MatrixXd elKe(3 * element->get_num_nodes(),3 * element->get_num_nodes());
    //        elKe.setZero();

    //        for(size_t g_id = 0; g_id < element->get_num_gaussian_points(); g_id++)
    //        {
    //            zxMaterialPoint::Ptr p_mat = element->get_material_point(g_id);

    //            zx_vega_ComputePartialPstress_to_DefGrad(m_material,p_mat->m_svd_U,p_mat->m_svd_diag,p_mat->m_svd_V,p_mat->m_dPsdF_diag);

    //            elKe += p_mat->m_dFdr_trans * p_mat->m_dPsdF_diag * p_mat->m_dFdr;
    //        }

    //    }

    tangentK *= 0.0;


    for(size_t el = 0; el < m_mesh->get_num_elements(); el++)
    {
        zxElement::Ptr element = m_mesh->get_element(el);

        Eigen::MatrixXd elKe(3 * element->get_num_nodes(),3 * element->get_num_nodes());
        elKe.setZero();

        for(size_t g_id = 0; g_id < element->get_num_gaussian_points(); g_id++)
        {
            zxMaterialPoint::Ptr p_mat = element->get_material_point(g_id);

            Eigen::MatrixXd& dFdx = p_mat->m_dFdr;//(9, 3 * neln);
            Eigen::MatrixXd& dPdF = p_mat->m_dPsdF_diag;//(9,9);
            Eigen::MatrixXd& dFdxT = p_mat->m_dFdr_trans;
            dPdF.setZero();
            double detJ0 = p_mat->m_detJac0;
            detJ0 *= element->get_gaussian_weight(g_id);
            zx_vega_ComputePartialPstress_to_DefGrad(m_material,p_mat->m_svd_U,p_mat->m_svd_diag,p_mat->m_svd_V,p_mat->m_dPsdF_diag);
            dPdF *= detJ0;
            elKe += dFdxT * dPdF.transpose() * dFdx;
        }

        for(size_t n_id0 = 0; n_id0 < element->get_num_nodes(); n_id0++)
            for(size_t n_id1 = 0; n_id1 < element->get_num_nodes(); n_id1++)
                for(int ii = 0; ii < 3; ii++)
                    for(int jj = 0; jj  < 3; jj++)
                    {
                        int r_id = 3 * element->get_node(n_id0)->m_id + ii;
                        int c_id = 3 * element->get_node(n_id1)->m_id + jj;

                        tangentK.coeffRef(r_id,c_id) += elKe(3 * n_id0 + ii, 3 * n_id1 + jj);
                    }
    }
}
