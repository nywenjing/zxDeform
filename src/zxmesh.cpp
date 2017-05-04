#include "zxmesh.h"

void zxSolidElement::init_material_points()
{
    for(size_t ni = 0; ni < get_num_gaussian_points(); ni++)
    {
        zxMaterialPoint::Ptr pt = get_material_point(ni);

        mat3d J;
        J.setZero();

        for (size_t nn=0; nn< get_num_nodes(); ++nn)
        {
            const double& Gri = get_first_derive_r(nn,ni);
            const double& Gsi = get_first_derive_s(nn,ni);
            const double& Gti = get_first_derive_t(nn,ni);

            const double& x = get_node(nn)->r0[0];
            const double& y = get_node(nn)->r0[1];
            const double& z = get_node(nn)->r0[2];

            J(0,0) += Gri*x; J(0,1) += Gsi*x; J(0,2) += Gti*x;
            J(1,0) += Gri*y; J(1,1) += Gsi*y; J(1,2) += Gti*y;
            J(2,0) += Gri*z; J(2,1) += Gsi*z; J(2,2) += Gti*z;
        }
        pt->m_invJac0 = J.inverse();
        pt->m_detJac0 = J.determinant();
    }

    for(size_t ni = 0; ni < get_num_gaussian_points(); ni++)
    {
        zxMaterialPoint::Ptr pt = get_material_point(ni);
        pt->m_dFdr.resize(9,3 * get_num_nodes());
        pt->m_dPsdF_diag.resize(9,9);


        Eigen::MatrixXd& dFdx = pt->m_dFdr;
        dFdx.setZero();
        mat3d& Ji0 = pt->m_invJac0;


        real* Gri = &get_first_derive_r(0,ni);
        real* Gsi = &get_first_derive_s(0,ni);
        real* Gti = &get_first_derive_t(0,ni);

        for(size_t nn = 0; nn < get_num_nodes(); nn++)
        {
            double Gxi = Ji0(0,0)*Gri[nn]+Ji0(1,0)*Gsi[nn]+Ji0(2,0)*Gti[nn];
            double Gyi = Ji0(0,1)*Gri[nn]+Ji0(1,1)*Gsi[nn]+Ji0(2,1)*Gti[nn];
            double Gzi = Ji0(0,2)*Gri[nn]+Ji0(1,2)*Gsi[nn]+Ji0(2,2)*Gti[nn];

            for(int ii = 0; ii < 3; ii++)
            {
                dFdx(ii * 3 + 0,3 * nn + ii) += Gxi;
                dFdx(ii * 3 + 1,3 * nn + ii) += Gyi;
                dFdx(ii * 3 + 2,3 * nn + ii) += Gzi;
            }
        }

        pt->m_dFdr_trans = dFdx.transpose();
    }
}

zxTetrahedralMesh::zxTetrahedralMesh(std::string filename)
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi T;

    std::string nodefile = filename + ".node";
    std::string elefile = filename + ".ele";


    std::ifstream ifile;
    ifile.open(nodefile);

    size_t nv;
    size_t dim;
    size_t flag0,flag1;
    size_t vid;

    ifile>>nv>>dim>>flag0>>flag1;

    V.resize(nv,3);
    for(size_t i = 0; i < nv; i++)
    {
        ifile>>vid >> V(i,0) >> V(i,1)>>V(i,2);
    }

    ifile.close();

    ifile.open(elefile);

    size_t nt;
    size_t el_num;
    size_t tag;
    size_t e_id;

    ifile>>nt>>el_num>>tag;

    T.resize(nt,4);

    for(size_t i = 0; i < nt; i++)
        ifile>>e_id>>T(i,0)>>T(i,1)>>T(i,2)>>T(i,3);


    ifile.close();


    m_nodes.resize(V.rows());
    m_elements.resize(T.rows());

    for(int n_id = 0; n_id < V.rows(); n_id++)
    {
        zxNode::Ptr node = zxNode::create();
        node->m_id = n_id;
        node->r0 = node->rp = node->rt = V.row(n_id);
        m_nodes[n_id] = node;
    }

    for(int e_id = 0; e_id < T.rows(); e_id++)
    {
        zxTetrahedron::Ptr tet = zxTetrahedron::create();
        tet->m_nodes.resize(T.cols());
        for(int j = 0; j < T.cols(); j++)
            tet->m_nodes[j] = get_node(T(e_id,j) - 1);

        tet->m_type = zxGaussianTrait::C3D4_1;
        tet->m_material_points.resize(tet->get_num_gaussian_points());
        for(int i = 0; i < tet->get_num_gaussian_points(); i++)
        {
            zxMaterialPoint::Ptr p_mat = zxMaterialPoint::create();
            p_mat->m_dPsdF_diag.resize(9,9);
            p_mat->m_dFdr.resize(9,12);
            p_mat->m_dFdr_trans.resize(12,9);
            tet->m_material_points[i] = p_mat;
        }

        tet->init_material_points();
        tet->m_id = e_id;
        m_elements[e_id] = tet;
    }
}

void zxTetrahedralMesh::buildSurface(Eigen::MatrixXd& V,Eigen::MatrixXi& W_id,Eigen::MatrixXd& W_val,Eigen::MatrixXi& F)
{
    class Face
    {
    public:
        size_t m_id[3];
        size_t m_id0[3];

    public:
        bool operator < (const Face& f2) const
        {
            for(int i = 0; i < 3; i++)
            {
                if(m_id[i] > f2.m_id[i])
                    return false;
                if(m_id[i] < f2.m_id[i])
                    return true;
            }

            return false;
        }

        bool operator == (const Face& f2) const
        {
            for(int i = 0; i < 3; i++)
                if(m_id[i] != f2.m_id[i])
                    return false;
            return true;
        }
    };

    class Tet
    {
    public:
        size_t m_id[4];
    };

    std::vector<Tet> allTets(get_num_elements());

    for(size_t i = 0; i < get_num_elements(); i++)
    {
        Tet& tet = allTets[i];
        for(int j = 0; j < 4; j++)
            tet.m_id[j] = get_node(i,j)->m_id;
    }

    std::vector<Face> allFaces(get_num_elements() * 4);

    for(size_t el = 0; el < get_num_elements(); el++)
    {
        Tet& tet = allTets[el];
        size_t face_indices[4][3] = {{0,2,1},{0,1,3},{0,3,2},{1,2,3}};

        for(int i = 0; i < 4; i++)
        {
            Face& face = allFaces[4 * el + i];
            for(int j = 0; j < 3; j++)
                face.m_id0[j] = face.m_id[j] = tet.m_id[face_indices[i][j]];

            std::sort(face.m_id,face.m_id + 3);
        }
    }

    std::sort(allFaces.begin(),allFaces.end());
    std::vector<Face> uniqFace = allFaces;
    std::vector<Face>::iterator last = std::unique_copy(allFaces.begin(),allFaces.end(),uniqFace.begin());

    uniqFace.resize(std::distance(uniqFace.begin(),last));

    std::vector<Face> singleFace;
    for(size_t i = 0; i < uniqFace.size(); i++)
    {
        Face& face = uniqFace[i];
        std::vector<Face>::iterator lb = std::lower_bound(allFaces.begin(),allFaces.end(),face);
        std::vector<Face>::iterator ub = std::upper_bound(allFaces.begin(),allFaces.end(),face);

        assert(lb != allFaces.end());

        if(std::distance(lb,ub) == 1)
            singleFace.push_back(face);

    }

    V.resize(get_num_nodes(),3);
    W_id.resize(get_num_nodes(),1);
    W_val.resize(get_num_nodes(),1);
    F.resize(singleFace.size(),3);
    for(size_t i = 0; i < get_num_nodes(); i++)
    {
        for(int j = 0; j < 3; j++)
            V(i,j) = get_node(i)->rt[j];
        W_id(i,0) = i;
        W_val(i,0) = 1.0;
    }

    for(size_t el = 0; el < singleFace.size(); el++)
    {
        for(int j = 0; j < 3; j++)
            F(el,j) = singleFace[el].m_id0[j];
    }
}
