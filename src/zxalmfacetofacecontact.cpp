#include "zxalmfacetofacecontact.h"
#include "igl/signed_distance.h"

void project_on_triangle(const vec3d &x0,
                         const vec3d &x1, const vec3d &x2, const vec3d &x3,double &s1, double &s2, double &s3)
{
    vec3d x13=x1-x3;
    double r00=(x13).norm()+1e-30;
    x13/=r00;
    vec3d x23=x2-x3;
    double r01=x23.dot(x13);
    x23-=r01*x13;
    double r11=(x23).norm()+1e-30;
    x23/=r11;
    vec3d x03=x0-x3;
    s2=x23.dot(x03)/r11;
    s1=(x13.dot(x03)-r01*s2)/r00;
    s3=1-s1-s2;
}

zxALMFaceToFaceContact::zxALMFaceToFaceContact(zxCollisionMesh::Ptr mastersurf,zxCollisionMesh::Ptr slavesurf,real eps)
    :zxALMContactInterface(eps)
{
    m_masterSurface = mastersurf;
    m_slaveSurface = slavesurf;

    m_gaussian_type = zxGaussianTrait::S3_3;

    m_almData.resize(m_slaveSurface->get_num_faces());

    m_num_gaussian = 0;
    for(size_t i = 0; i < m_almData.size(); i++)
    {
        m_almData[i].resize(zxGaussianFactory::Singleton().get_num_gaussian_points(m_gaussian_type));

        for(size_t j = 0; j < m_almData[i].size(); j++)
        {
            ALMData& data = m_almData[i][j];
            data.m_id = m_num_gaussian++;
        }
    }
}


void   zxALMFaceToFaceContact::update_contact()
{
    Eigen::MatrixXd P,V;
    Eigen::MatrixXi F;
    Eigen::VectorXd S;
    Eigen::VectorXi I;
    Eigen::MatrixXd C,N;

    m_slaveSurface->update_position();
    m_masterSurface->update_position();

    P.resize(m_num_gaussian,3);

    for(size_t el = 0; el < m_slaveSurface->get_num_faces(); el++)
        for(size_t g_id = 0; g_id < zxGaussianFactory::Singleton().get_num_gaussian_points(m_gaussian_type); g_id++)
        {
            vec3d x = vec3d::Zero();

            for(size_t n_id = 0; n_id < 3; n_id++)
                x += zxGaussianFactory::Singleton().get_shape_fun(m_gaussian_type,n_id,g_id) * m_slaveSurface->get_face(el)->v[n_id]->x;


            for(size_t k = 0; k < 3; k++)
                P(m_almData[el][g_id].m_id,k) = x[k];
            m_almData[el][g_id].m_x = x;
        }

    m_masterSurface->get_Matrix_Format(V,F);

    igl::signed_distance(P,V,F,igl::SignedDistanceType::SIGNED_DISTANCE_TYPE_UNSIGNED,
                         S,I,C,N);

    for(size_t el = 0; el < m_slaveSurface->get_num_faces(); el++)
        for(size_t g_id = 0; g_id < zxGaussianFactory::Singleton().get_num_gaussian_points(m_gaussian_type); g_id++)
        {
            ALMData& data = m_almData[el][g_id];
            zxCollisionMesh::Face::Ptr face = m_masterSurface->get_face(I[m_almData[el][g_id].m_id]);
            double s0,s1,s2;

            project_on_triangle(m_almData[el][g_id].m_x,face->v[0]->x,face->v[1]->x, face->v[2]->x,s0,s1,s2);

            vec3d q = s0 * face->v[0]->x + s1 * face->v[1]->x + s2 * face->v[2]->x;

            data.m_master_element = nullptr;

            if(s0 >= zxEPSILON && s1 >= zxEPSILON && s2 >= zxEPSILON )
            {
                data.m_master_element = face.get();
                data.m_normal = face->compute_normal();
                data.m_gap = data.m_normal.dot(q - data.m_x);
                data.m_r = s1;
                data.m_s = s2;
            }


        }

}

void   zxALMFaceToFaceContact::compute_contact_force(Eigen::VectorXd& R)
{
    vec3d mTotalActionForce(0.0,0.0,0.0);
    for(size_t el = 0; el < m_slaveSurface->get_num_faces(); el++)
    {
        zxCollisionMesh::Face* se = m_slaveSurface->get_face(el).get();
        real area = se->compute_area0();
        for(size_t g_id = 0; g_id < zxGaussianFactory::Singleton().get_num_gaussian_points(m_gaussian_type); g_id++)
        {
            ALMData& data = m_almData[el][g_id];


            if(data.m_master_element == nullptr)
                continue;

            real tn = data.m_lambda + m_eps * data.m_gap;
            real gw = zxGaussianFactory::Singleton().get_gaussian_weight(m_gaussian_type,g_id);


            if(tn < 0) tn = 0;
            if(tn  == 0) continue;

            vec3d traction = data.m_normal * tn * area * 2.0 * gw;
            zxCollisionMesh::Face* me = data.m_master_element;


            mTotalActionForce += traction;
            real Hm3[3],Hs3[3];

            zxGaussianFactory::Singleton().get_shape_fun(m_gaussian_type,Hm3,data.m_r,data.m_s,0.0);
            for(size_t nid = 0; nid < 3; nid++)
                Hs3[nid] = zxGaussianFactory::Singleton().get_shape_fun(m_gaussian_type,nid,g_id);


            for(size_t n_id = 0; n_id < 3; n_id++)
            {
                real Hm = Hm3[n_id];

                for(size_t i = 0; i <me->v[n_id]->m_nodes.size(); i++)
                    for(size_t j = 0; j < 3; j++)
                    {
                        int eid = me->v[n_id]->m_nodes[i]->m_dof_id[j];
                        if(eid >= 0)
                            R[eid] -= traction[j] * Hm * me->v[n_id]->m_weights[i];
                    }
            }

            for(size_t n_id = 0; n_id < 3; n_id++)
            {
                real Hs = Hs3[n_id];

                for(size_t i = 0; i <se->v[n_id]->m_nodes.size(); i++)
                    for(size_t j = 0; j < 3; j++)
                    {
                        int eid = se->v[n_id]->m_nodes[i]->m_dof_id[j];
                        if(eid >= 0)
                            R[eid] += traction[j] * Hs * se->v[n_id]->m_weights[i];
                    }
            }

        }
    }

    //std::cout<<"contact force\n"<<mTotalActionForce<<std::endl;

}

std::list<Eigen::Triplet<real>>   zxALMFaceToFaceContact::get_stiffness_triplets()
{
    std::list<Eigen::Triplet<real>> stiffEntries;

    for(size_t el = 0; el < m_slaveSurface->get_num_faces(); el++)
        for(size_t g_id = 0; g_id < zxGaussianFactory::Singleton().get_num_gaussian_points(m_gaussian_type); g_id++)
        {
            ALMData& data = m_almData[el][g_id];

            if(data.m_master_element == nullptr)
                continue;
            zxCollisionMesh::Face* me = data.m_master_element;
            zxCollisionMesh::Face* se = m_slaveSurface->get_face(el).get();
            for(size_t i = 0; i < 3; i++)
                for(size_t ii = 0 ; ii < se->v[i]->m_nodes.size(); ii++)
                    for(size_t j = 0; j < 3; j++)
                        for(size_t jj = 0; jj < me->v[j]->m_nodes.size();jj++)
                        {
                            zxNode* mnode = me->v[j]->m_nodes[jj].get();
                            zxNode* snode = se->v[i]->m_nodes[ii].get();

                            for(size_t l = 0; l < 3; l++)
                                for(size_t m = 0; m < 3; m++)
                                {
                                    int dof0 = mnode->m_dof_id[l];
                                    int dof1 = snode->m_dof_id[m];

                                    if(dof0 < -1)
                                        dof0 = -dof0 - 2;
                                    if(dof1 < -1)
                                        dof1 = -dof1 - 2;
                                    if(dof0 >= 0 && dof1 >= 0)
                                        stiffEntries.push_back(Eigen::Triplet<real>(dof0,dof1,0.0));
                                }
                        }
        }

    return stiffEntries;
}

void   zxALMFaceToFaceContact::compute_stiffness(Eigen::SparseMatrix<real>& tangentStiff)
{
    std::list<Eigen::Triplet<real>> stiffEntries;

    std::vector<int>    dofvec;
    std::vector<real>   weightvec;

    dofvec.reserve(3 * (6 + 6));
    weightvec.reserve(3 * (6 + 6));

    for(size_t el = 0; el < m_slaveSurface->get_num_faces(); el++)
    {
        zxCollisionMesh::Face* se = m_slaveSurface->get_face(el).get();
        real area = se->compute_area0();

        for(size_t g_id = 0; g_id < zxGaussianFactory::Singleton().get_num_gaussian_points(m_gaussian_type); g_id++)
        {
            dofvec.clear();
            weightvec.clear();

            ALMData& data = m_almData[el][g_id];

            if(data.m_master_element == nullptr)
                continue;

            real tn = data.m_lambda + m_eps * data.m_gap;
            real gw = zxGaussianFactory::Singleton().get_gaussian_weight(m_gaussian_type,g_id);

            if(tn < 0)          tn = 0;
            if(tn == 0)  continue;
            vec3d normal = data.m_normal;
            real Hm3[3];
            real Hs3[3];

            zxGaussianFactory::Singleton().get_shape_fun(m_gaussian_type,Hm3,data.m_r,data.m_s,0.0);
            for(size_t nid = 0; nid < 3; nid++)
                Hs3[nid] = zxGaussianFactory::Singleton().get_shape_fun(m_gaussian_type,nid,g_id);


            zxCollisionMesh::Face* me = data.m_master_element;

            for(size_t i = 0; i < 3; i++)
                for(size_t j = 0; j < se->v[i]->m_nodes.size(); j++)
                {
                    zxNode* node = se->v[i]->m_nodes[j].get();

                    for(size_t k = 0; k < 3; k++)
                    {
                        dofvec.push_back(node->m_dof_id[k]);
                        weightvec.push_back(Hs3[i] * se->v[i]->m_weights[j] * normal[k]);
                    }
                }

            for(size_t i = 0; i < 3; i++)
                for(size_t j = 0; j < me->v[i]->m_nodes.size(); j++)
                {
                    zxNode* node = me->v[i]->m_nodes[j].get();

                    for(size_t k = 0; k < 3; k++)
                    {
                        dofvec.push_back(node->m_dof_id[k]);
                        weightvec.push_back(-Hm3[i] * me->v[i]->m_weights[j] * normal[k]);
                    }
                }

            for(size_t i = 0; i < dofvec.size(); i++)
                for(size_t j = 0; j < dofvec.size(); j++)
                    if(dofvec[i] >= 0 && dofvec[j] >= 0)
                    {
                        real kij = weightvec[i] * weightvec[j];
                        kij *= area * 2.0 * gw * m_eps;

                        int dof0 = dofvec[i];
                        int dof1 = dofvec[j];

                        if(dof0 < -1)
                            dof0 = -dof0 - 2;
                        if(dof1 < -1)
                            dof1 = -dof1 - 2;
                        if(dof0 >= 0 && dof1 >= 0)


                        stiffEntries.push_back(Eigen::Triplet<real>(dof0,dof1,gw * area * m_eps * 2 * weightvec[i] * weightvec[j]));
                    }

        }
    }

    tangentStiff.setFromTriplets(stiffEntries.begin(),stiffEntries.end());
}

bool   zxALMFaceToFaceContact::alm_augment()
{
    real lnorm0 = 0;
    real lnorm1 = 0;

    for(size_t el = 0; el < m_almData.size(); el++)
        for(size_t g_id = 0; g_id < m_almData[el].size(); g_id++)
        {
            ALMData& data = m_almData[el][g_id];
            lnorm0 += data.m_lambda * data.m_lambda;
        }

    for(size_t el = 0; el < m_almData.size(); el++)
        for(size_t g_id = 0; g_id < m_almData[el].size(); g_id++)
        {
            ALMData& data = m_almData[el][g_id];
            real lambda = data.m_lambda + m_eps * data.m_gap;

            if(lambda < 0) lambda = 0;
            lnorm1 += lambda * lambda;
        }

    lnorm0 = std::sqrt(lnorm0);
    lnorm1 = std::sqrt(lnorm1);

    real cratio = 0;
    if(lnorm1 > zxEPSILON)
        cratio = std::abs(lnorm1 - lnorm0) / lnorm1;

    bool converge = cratio < m_alm_tol;

    if(!converge)
    {
        for(size_t el = 0; el < m_almData.size(); el++)
            for(size_t g_id = 0; g_id < m_almData[el].size(); g_id++)
            {
                ALMData& data = m_almData[el][g_id];
                real lambda = data.m_lambda + m_eps * data.m_gap;

                if(lambda < 0) lambda = 0;
                data.m_lambda = lambda;
            }
    }
    return converge;
}
