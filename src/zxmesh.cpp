#include "zxmesh.h"

void zxSolidElement::init_material_points()
{
    m_material_points.resize(get_num_gaussian_points());
    for(int i = 0; i < get_num_gaussian_points(); i++)
    {
        m_material_points[i] = zxMaterialPoint::create();
    }

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

        assert(pt->m_detJac0 > 0);
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

    assert(ifile.is_open());

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


        tet->init_material_points();
        tet->m_id = e_id;
        m_elements[e_id] = tet;
    }
}

void zxTetrahedralMesh::convertToC3D10()
{
    assert(get_element(0)->m_nodes.size() != 10);

    class Edge
    {
    public:
        size_t m_id[2];
        int    m_gid;

        zxNode::Ptr m_node;

        bool operator < (const Edge& e) const
        {
            for(int i = 0; i < 2; i++)
            {
                if(m_id[i] > e.m_id[i])
                    return false;
                if(m_id[i] < e.m_id[i])
                    return true;
            }

            return false;

        }

        bool operator == (const Edge& e) const
        {
            return (m_id[0] == e.m_id[0]) && (m_id[1] == e.m_id[1]);
        }
    };

    int elEdge[6][2] = {{0,1},{1,2},{0,2},{0,3},{1,3},{2,3}};

    std::vector<Edge> allEdges;
    allEdges.reserve(6 * get_num_elements());

    for(size_t el = 0; el < get_num_elements(); el++)
    {
        for(int i = 0; i < 6; i++)
        {
            Edge e;
            int id0 = get_node(el,elEdge[i][0])->m_id;
            int id1 = get_node(el,elEdge[i][1])->m_id;

            e.m_id[0] = std::min(id0,id1);
            e.m_id[1] = std::max(id0,id1);

            allEdges.push_back(e);
        }
    }

    std::sort(allEdges.begin(),allEdges.end());

    std::vector<Edge> uniqEdge = allEdges;
    std::vector<Edge>::iterator last = std::unique_copy(allEdges.begin(),allEdges.end(),uniqEdge.begin());
    uniqEdge.resize(std::distance(uniqEdge.begin(),last));

    m_nodes.reserve(m_nodes.size() + uniqEdge.size());
    for(size_t i = 0; i < uniqEdge.size(); i++)
    {
        Edge& e = uniqEdge[i];
        e.m_gid = i;
        zxNode::Ptr n0 = get_node(e.m_id[0]);
        zxNode::Ptr n1 = get_node(e.m_id[1]);

        zxNode::Ptr n_mid = zxNode::create();
        n_mid->r0 = 0.5 * (n0->r0 + n1->r0);
        n_mid->rp = 0.5 * (n0->rp + n1->rp);
        n_mid->rt = 0.5 * (n0->rt + n1->rt);

        n_mid->m_id = i + get_num_nodes();

        e.m_node = n_mid;


    }

    for(size_t i = 0; i < uniqEdge.size(); i++)
    {
        Edge& e = uniqEdge[i];
        m_nodes.push_back(e.m_node);
    }

    for(size_t el = 0; el < get_num_elements(); el++)
    {
        for(int i = 0; i < 6; i++)
        {
            Edge e;
            int id0 = get_node(el,elEdge[i][0])->m_id;
            int id1 = get_node(el,elEdge[i][1])->m_id;

            e.m_id[0] = std::min(id0,id1);
            e.m_id[1] = std::max(id0,id1);

            std::vector<Edge>::iterator lb = std::lower_bound(uniqEdge.begin(),uniqEdge.end(),e);

            get_element(el)->m_nodes.push_back((*lb).m_node);
        }

        get_element(el)->m_type = zxGaussianTrait::C3D10_4;

        zxTetrahedron* tet = static_cast<zxTetrahedron*>(get_element(el).get());
        tet->init_material_points();
    }


}

void zxTetrahedralMesh::buildSurface(Eigen::MatrixXd& V,Eigen::MatrixXi& W_id,Eigen::MatrixXd& W_val,Eigen::MatrixXi& F,int tsLevel)
{
    class Face
    {
    public:
        size_t m_id[3];
        size_t m_id0[3];
        size_t m_tid;
        size_t m_tfid;

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

    for(size_t el = 0; el < get_num_elements(); el++)
    {
        Tet& tet = allTets[el];
        size_t face_indices[4][3] = {{0,2,1},{0,1,3},{0,3,2},{1,2,3}};

        for(int i = 0; i < 4; i++)
        {
            Face face;
            for(int j = 0; j < 3; j++)
                face.m_id0[j] = face.m_id[j] = tet.m_id[face_indices[i][j]];

            std::sort(face.m_id,face.m_id + 3);

            std::vector<Face>::iterator lb = std::lower_bound(singleFace.begin(),singleFace.end(),face);


            if(lb != singleFace.end() && (*lb) == face)
            {
                (*lb).m_tid = el;
                (*lb).m_tfid = i;
            }
        }
    }

    if(get_element(0)->m_type == zxGaussianTrait::C3D4_1)
    {
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
    else if(get_element(0)->m_type == zxGaussianTrait::C3D10_4)
    {
        class LevelFace
        {
            ZX_MAKE_SHARED_MACO_DEFAULT_CREAT(LevelFace)
                    public:
                vec2d m_rs[3];
            size_t m_level;
            LevelFace::Ptr m_child[4];

            std::vector<LevelFace::Ptr> m_allLeaf;
        public:
            void setLevel(int le,vec2d rs0,vec2d rs1,vec2d rs2)
            {
                m_level = le;
                m_rs[0] = rs0;
                m_rs[1] = rs1;
                m_rs[2] = rs2;

                vec2d rs3 = 0.5 * (rs0 + rs1);
                vec2d rs4 = 0.5 * (rs1 + rs2);
                vec2d rs5 = 0.5 * (rs0 + rs2);

                if(le != 0)
                {
                    m_child[0] = LevelFace::create();
                    m_child[0]->setLevel(le-1,rs0,rs3,rs5);

                    m_child[1] = LevelFace::create();
                    m_child[1]->setLevel(le-1,rs3,rs1,rs4);

                    m_child[2] = LevelFace::create();
                    m_child[2]->setLevel(le-1,rs5,rs4,rs2);

                    m_child[3] = LevelFace::create();
                    m_child[3]->setLevel(le-1,rs3,rs4,rs5);

                }
            }

            void collectLeaf()
            {
                for(int i = 0; i < 4; i++)
                {
                    if(m_child[i]->m_level > 0)
                    {
                        m_child[i]->collectLeaf();
                        m_allLeaf.insert(m_allLeaf.end(),m_child[i]->m_allLeaf.begin(),m_child[i]->m_allLeaf.end());
                    }
                    else
                        m_allLeaf.push_back(m_child[i]);
                }
            }

        };

        LevelFace::Ptr lface = LevelFace::create();
        lface->setLevel(tsLevel + 1,vec2d(0,0),vec2d(1,0),vec2d(0,1));
        lface->collectLeaf();

        class tVert
        {
        public:
            tVert()
            {
                this->m_w.resize(6);
                this->m_node.resize(6);
            }

            bool operator < (const tVert& vert) const
            {
                if(*this == vert)
                    return false;
                for(int i = 0; i < 3; i++)
                {
                    if(x[i] > vert.x[i])
                        return false;

                    if(x[i] < vert.x[i])
                        return true;
                }

                return false;
            }

            bool operator == (const tVert& vert) const
            {
                return (x - vert.x).norm() < zxEPSILON;
            }
        public:
            int   m_id;
            vec3d x;
            std::vector<zxNode::Ptr> m_node;
            std::vector<real>    m_w;
        };

        std::list<tVert> tVerts;


        for(size_t el = 0; el < singleFace.size(); el++)
        {
            Face& face = singleFace[el];

            size_t faceIdx[4][6] = {
                {0,2,1,6,5,4},
                {0,1,3,4,8,7},
                {0,3,2,7,9,6},
                {1,2,3,5,9,8}
            };

            std::vector<zxNode::Ptr> fnode(6);
            for(int i = 0; i < 6; i++)
                fnode[i] = get_node(face.m_tid,faceIdx[face.m_tfid][i]);


            for(size_t se = 0 ; se < lface->m_allLeaf.size(); se++)
            {
                LevelFace::Ptr sface = lface->m_allLeaf[se];


                for(size_t i = 0; i < 3; i++)
                {
                    tVert tv;
                    zxGaussianFactory::Singleton().get_shape_fun(zxGaussianTrait::S6_3,tv.m_w.data(),sface->m_rs[i][0],sface->m_rs[i][1],0.0);
                    tv.m_node = fnode;

                    tv.x = vec3d::Zero();
                    for(int j = 0; j < 6; j++)
                        tv.x += tv.m_w[j] * tv.m_node[j]->rt;

                    tVerts.push_back(tv);
                }

            }


        }

        V.resize(tVerts.size(),3);
        W_id.resize(tVerts.size(),6);
        W_val.resize(tVerts.size(),6);
        F.resize(tVerts.size()/3,3);

        int v_id = 0;
        int f_id = 0;

        std::vector<tVert> vectVerts(tVerts.begin(),tVerts.end());
        for(int i = 0; i < vectVerts.size(); i++)
            vectVerts[i].m_id = i;

        std::vector<tVert> sort_verts = vectVerts;
        std::sort(sort_verts.begin(),sort_verts.end());
        std::vector<tVert> uniq_vert = sort_verts;
        std::vector<tVert>::iterator last = std::unique_copy(sort_verts.begin(),sort_verts.end(),uniq_vert.begin());
        uniq_vert.resize(std::distance(uniq_vert.begin(),last));

        V.resize(uniq_vert.size(),3);
        W_id.resize(uniq_vert.size(),6);
        W_val.resize(uniq_vert.size(),6);

        for(size_t i = 0; i < uniq_vert.size(); i++)
        {
            tVert& tv = uniq_vert[i];
            for(int j = 0; j < 3; j++)
                V(i,j) = uniq_vert[i].x[j];

            for(int j = 0; j < 6; j++)
            {
                W_id(i,j) = tv.m_node[j]->m_id;
                W_val(i,j) = tv.m_w[j];
            }

            tv.m_id = i;
        }

        for(int el = 0; el < F.rows(); el++)
            for(int i = 0; i < 3; i++)
            {
                std::vector<tVert>::iterator lb = std::lower_bound(uniq_vert.begin(),uniq_vert.end(),vectVerts[3 * el + i]);

                F(el,i) = (*lb).m_id;
            }



//        for(std::list<tVert>::iterator it = tVerts.begin();it != tVerts.end();)
//        {
//            for(int vi = 0; vi < 3; vi++)
//            {
//                tVert& tv = *it++;
//                for(int i = 0; i < 3; i++)
//                    V(v_id,i) = tv.x[i];
//                for(int i = 0; i < 6; i++)
//                {
//                    W_id(v_id,i) = tv.m_node[i]->m_id;
//                    W_val(v_id,i) = tv.m_w[i];
//                }

//                F(f_id,vi) = v_id;

//                v_id++;
//            }

//            f_id++;
//        }

        //        V.resize(singleFace.size() * 3, 3);
        //        W_id.resize(singleFace.size() * 3,3);
        //        W_val.resize(singleFace.size() * 3,3);
        //        F.resize(singleFace.size(),3);

        //        W_val.setZero();

        //        for(size_t el = 0; el < singleFace.size(); el++)
        //        {
        //            Face& face = singleFace[el];

        //            size_t faceIdx[4][6] = {
        //                {0,2,1,6,5,4},
        //                {0,1,3,4,8,7},
        //                {0,3,2,7,9,6},
        //                {1,2,3,5,9,8}
        //            };

        //            std::vector<zxNode::Ptr> fnode(6);
        //            for(int i = 0; i < 6; i++)
        //                fnode[i] = get_node(face.m_tid,faceIdx[face.m_tfid][i]);

        //            int fid[6];
        //            for(int i = 0; i < 6; i++)
        //                fid[i] = fnode[i]->m_id;
        //            int tid[10];
        //            for(int i = 0; i < 10; i++)
        //                tid[i] = get_node(face.m_tid,i)->m_id;

        //            for(int vi = 0; vi < 3; vi++)
        //            {
        //                F(el,vi) = 3 * el + vi;

        //                for(int i = 0; i < 3; i++)
        //                {
        //                    W_id(3 * el + vi,i) = fnode[i]->m_id;//face.m_id0[i];
        //                    W_val(3 * el + vi,vi) = 1.0;
        //                }
        //            }

        //        }



    }


}
