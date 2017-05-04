#ifndef ZXCOLLISIONMESH_H
#define ZXCOLLISIONMESH_H

#include "zxmesh.h"
#include "zxbvh.h"

class zxCollisionMesh
{
    ZX_MAKE_SHARED_MACO(zxCollisionMesh)

public:
     static Ptr create(std::string filename){ return Ptr(new zxCollisionMesh(filename));}
    static Ptr create(zxSolidMesh::Ptr mesh){ return Ptr(new zxCollisionMesh(mesh));}
     zxCollisionMesh(std::string filename);
     zxCollisionMesh(zxSolidMesh::Ptr mesh);

public:
    class Face;
    class Vert
    {
        ZX_MAKE_SHARED_MACO_DEFAULT_CREAT(Vert)
    public:
        void update_position()
        {
            x = xp = x0 = vec3d::Zero();
            for(size_t i = 0; i < m_nodes.size(); i++)
            {
                zxNode::Ptr node = m_nodes[i];

                x += node->rt * m_weights[i];
                xp += node->rp * m_weights[i];
                x0 += node->r0 * m_weights[i];
            }

        }

    public:
        std::vector<zxNode::Ptr> m_nodes;
        std::vector<real>        m_weights;
    public:
        vec3d x,xp,x0;
        int   m_id;
        Face* m_test_face;
    };

    class Edge
    {
        ZX_MAKE_SHARED_MACO_DEFAULT_CREAT(Edge)

    public:
        Vert::Ptr   v[2];
        Face* m_test_face;
    };

    class Face : public zxAABBData
    {
        ZX_MAKE_SHARED_MACO_DEFAULT_CREAT(Face)


        public:
        virtual void update_aabb(bool ccd)
        {
            m_aabb->reset();

            for(size_t i = 0; i < 3; i++)
                m_aabb->combine(v[i]->x);

            if(ccd)
            {
                for(size_t i = 0; i < 3; i++)
                    m_aabb->combine(v[i]->xp);
            }

        }

        public:
            Vert::Ptr   v[3];
            Edge::Ptr   e[3];
    };

public:
    zxCollisionMesh();
    virtual void build_edges();
public:
    size_t      get_num_verts() {return m_verts.size();}
    size_t      get_num_edges() {return m_edges.size();}
    size_t      get_num_faces() {return m_faces.size();}

    Vert::Ptr   get_vert(int id){return m_verts.at(id);}
    Edge::Ptr   get_edge(int id){return m_edges.at(id);}
    Face::Ptr   get_face(int id){return m_faces.at(id);}

    std::vector<Face::Ptr>& get_faces(){return m_faces;}

public:
    void        get_Matrix_Format(Eigen::MatrixXd& V,Eigen::MatrixXi& F);

public:
    void        update_aabb(bool ccd);
    void        update_position();



protected:

    std::vector<Vert::Ptr>  m_verts;
    std::vector<Edge::Ptr>  m_edges;
    std::vector<Face::Ptr>  m_faces;


};

#endif // ZXCOLLISIONMESH_H