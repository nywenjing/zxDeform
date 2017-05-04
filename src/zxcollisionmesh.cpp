#include "zxcollisionmesh.h"
#include "igl/read_triangle_mesh.h"
zxCollisionMesh::zxCollisionMesh()
{

}

zxCollisionMesh::zxCollisionMesh(std::string filename)
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh(filename,V,F);

    m_verts.resize(V.rows());
    m_faces.resize(F.rows());

    for(int i = 0; i < m_verts.size(); i++)
    {
        Vert::Ptr v = Vert::create();
        v->x = v->x0 = v->xp = V.row(i);
        v->m_id = i;

        m_verts[i] = v;
    }

    for(int el = 0; el < m_faces.size(); el++)
    {
        Face::Ptr face = Face::create();

        for(int j = 0; j < 3; j++)
            face->v[j] = get_vert(F(el,j));

        m_faces[el] = face;
    }

    build_edges();
}

zxCollisionMesh::zxCollisionMesh(zxSolidMesh::Ptr mesh)
{
    Eigen::MatrixXd V,W_val;
    Eigen::MatrixXi F,W_id;
    mesh->buildSurface(V,W_id,W_val,F);

    m_verts.resize(V.rows());
    m_faces.resize(F.rows());

    for(int i = 0; i < m_verts.size(); i++)
    {
        Vert::Ptr v = Vert::create();
        v->x = v->x0 = v->xp = V.row(i);
        v->m_id = i;

        m_verts[i] = v;

        for(int j = 0; j < W_id.cols(); j++)
        {
            v->m_nodes.push_back(mesh->get_node(W_id(i,j)));
            v->m_weights.push_back(W_val(i,j));
        }
    }

    for(int el = 0; el < m_faces.size(); el++)
    {
        Face::Ptr face = Face::create();

        for(int j = 0; j < 3; j++)
            face->v[j] = get_vert(F(el,j));

        m_faces[el] = face;
    }

    build_edges();

}

void zxCollisionMesh::build_edges()
{
    class tEdge
    {
    public:
        int id[2];
        bool operator < (const tEdge& e) const
        {
            for(int i = 0; i < 2; i++)
            {
                if(id[i] < e.id[i])
                    return true;
                if(id[i] > e.id[i])
                    return false;
            }

            return false;
        }

        bool operator == (const tEdge& e) const
        {
            return (id[0] == e.id[0]) && (id[1] == e.id[1]);
        }
    };

    std::vector<tEdge> all_edges(get_num_faces() * 3);

    for(size_t el = 0; el < get_num_faces(); el++)
    {
        Face::Ptr face = get_face(el);
        for(int i = 0; i < 3; i++)
        {
            int id0 = face->v[i]->m_id;
            int id1 = face->v[(i+1)%3]->m_id;

            tEdge& edge = all_edges[3 * el + i];
            edge.id[0] = std::min(id0,id1);
            edge.id[1] = std::max(id0,id1);
        }
    }

    std::sort(all_edges.begin(),all_edges.end());
    std::vector<tEdge> single_edges = all_edges;
    std::vector<tEdge>::iterator last = std::unique_copy(all_edges.begin(),all_edges.end(),single_edges.begin());
    single_edges.resize(std::distance(single_edges.begin(),last));

    m_edges.resize(single_edges.size());

    for(int i = 0; i < m_edges.size(); i++)
    {
        Edge::Ptr edge = Edge::create();
        tEdge& tedge = single_edges[i];
        edge->v[0] = get_vert(tedge.id[0]);
        edge->v[1] = get_vert(tedge.id[1]);
        m_edges[i] = edge;
    }

    for(size_t el = 0; el < get_num_faces(); el++)
    {
        Face::Ptr face = get_face(el);
        for(int i = 0; i < 3; i++)
        {
            int id0 = face->v[i]->m_id;
            int id1 = face->v[(i+1)%3]->m_id;

            tEdge edge;
            edge.id[0] = std::min(id0,id1);
            edge.id[1] = std::max(id0,id1);

            std::vector<tEdge>::iterator iter = std::lower_bound(single_edges.begin(),single_edges.end(),edge);

            assert(iter != single_edges.end());

            size_t e_id = std::distance(single_edges.begin(),iter);
            face->e[i] = get_edge(e_id);
        }
    }


    for(size_t el = 0; el < get_num_faces(); el++)
    {
        Face::Ptr face = get_face(el);

        for(size_t i = 0; i < 3; i++)
        {
            face->v[i]->m_test_face = face.get();
            face->e[i]->m_test_face = face.get();
        }
    }

}

void zxCollisionMesh::get_Matrix_Format(Eigen::MatrixXd& V,Eigen::MatrixXi& F)
{
    V.resize(get_num_verts(),3);
    F.resize(get_num_faces(),3);

    for(size_t i = 0; i < get_num_verts(); i++)
    {
        Vert::ConstPtr vert = get_vert(i);
        for(size_t j = 0; j < 3; j++)
            V(i,j) = vert->x[j];
    }

    for(size_t el = 0; el < get_num_faces(); el++)
        for(size_t j = 0; j < 3; j++)
            F(el,j) = get_face(el)->v[j]->m_id;

}

void zxCollisionMesh::update_aabb(bool ccd)
{
    for(size_t i = 0; i < get_num_faces(); i++)
    {
        Face::Ptr face = get_face(i);
        face->update_aabb(ccd);
    }
}

void zxCollisionMesh::update_position()
{
    for(size_t i = 0; i < get_num_verts(); i++)
    {
        Vert::Ptr vert = get_vert(i);
        vert->update_position();
    }
}
