#include "zxcollisionalgorithm.h"
#include "zxmath.h"
#include <memory>
zxCollisionAlgorithm::zxCollisionAlgorithm()
{

}

void zxFaceFace_Collider::process_collision(zxBVHNode* node0, zxBVHNode* node1)
{

    zxCollisionMesh::Face::Ptr face0 = std::static_pointer_cast<zxCollisionMesh::Face> (node0->get_data());
    zxCollisionMesh::Face::Ptr face1 = std::static_pointer_cast<zxCollisionMesh::Face> (node1->get_data());

    assert(face0 != nullptr);
    assert(face1 != nullptr);

    for(size_t i = 0; i < 3; i++)
        for(size_t j = 0; j < 3; j++)
        {
            if(face0->v[i] == face1->v[j])
                return;
        }

    for(size_t i = 0; i < 3; i++)
    {
        zxCollisionMesh::Vert::Ptr vert = face0->v[i];
        if(vert->m_test_face == face0.get())
            do_VF_Test(vert,face1);
    }

    for(size_t i = 0; i < 3; i++)
    {
        zxCollisionMesh::Vert::Ptr vert = face1->v[i];

        if(vert->m_test_face == face1.get())
            do_VF_Test(vert,face0);
    }

    for(size_t i = 0; i < 3; i++)
        for(size_t j = 0; j < 3; j++)
        {
            zxCollisionMesh::Edge::Ptr e0 = face0->e[i];
            zxCollisionMesh::Edge::Ptr e1 = face1->e[j];

            if((e0->m_test_face == face0.get()) && (e1->m_test_face == face1.get()))
                do_EE_Test(e0,e1);
        }
}

void zxFaceFace_Proximity::do_VF_Test(zxCollisionMesh::Vert::Ptr vert, zxCollisionMesh::Face::Ptr face)
{
    real s1,s2,s3;
    vec3d normal;
    real distance;

    zx_check_point_triangle_proximity(
                vert->x,
                face->v[0]->x,face->v[1]->x,face->v[2]->x,
            distance,
            s1,s2,s3,
            normal
            );


    bool checked = distance < m_margin && s1 >= -zxEPSILON && s2 >= -zxEPSILON && s3 >= -zxEPSILON;
    if(checked)
    {

        //std::cout<<"test1: "<<t<<" "<<s1<<" "<<s2<<" "<<s3<<std::endl;
        zxContactPoint::Ptr contact = zxContactPoint::create();
        contact->m_normal = normal;
        contact->m_verts.resize(4);
        contact->m_weights.resize(4);
        contact->m_verts[0] = vert;
        contact->m_verts[1] = face->v[0];
        contact->m_verts[2] = face->v[1];
        contact->m_verts[3] = face->v[2];

        contact->m_weights[0] = 1.0;
        contact->m_weights[1] = -s1;
        contact->m_weights[2] = -s2;
        contact->m_weights[3] = -s3;
        contact->m_thickness = m_margin * 0.9;
        m_contacts.push_back(contact);
    }



}

void zxFaceFace_Proximity::do_EE_Test(zxCollisionMesh::Edge::Ptr edge0, zxCollisionMesh::Edge::Ptr edge1)
{
    double EPSILON = 1e-6;

    real s0,s1,s2,s3;
    vec3d normal;
    real distance;

    zx_check_edge_edge_proximity(
                edge0->v[0]->x,edge0->v[1]->x,
            edge1->v[0]->x,edge1->v[1]->x,
            distance,
            s0,s2,
            normal);

    bool checked = false;
    checked = distance < m_margin && s0 > 0 && s0 < 1 && s2 > 0 && s2 < 1;

    s1 = 1 - s0;
    s3 = 1 - s2;

    if(checked)
    {

        zxContactPoint::Ptr contact = zxContactPoint::create();
        contact->m_normal = normal;
        contact->m_verts.resize(4);
        contact->m_weights.resize(4);
        contact->m_verts[0] = edge0->v[0];
        contact->m_verts[1] = edge0->v[1];
        contact->m_verts[2] = edge1->v[0];
        contact->m_verts[3] = edge1->v[1];

        contact->m_weights[0] = s0;
        contact->m_weights[1] = s1;
        contact->m_weights[2] = -s2;
        contact->m_weights[3] = -s3;
        contact->m_thickness = m_margin * 0.9;
        m_contacts.push_back(contact);
    }

}

void zxFaceFace_CCD::do_VF_Test(zxCollisionMesh::Vert::Ptr vert, zxCollisionMesh::Face::Ptr face)
{

    real s1,s2,s3;
    vec3d normal;
    real t;
    real collision_epsilon = zxEPSILON;

    bool checked = zx_check_point_triangle_collision(vert->xp, face->v[0]->xp, face->v[1]->xp, face->v[2]->xp,
            vert->x, face->v[0]->x, face->v[1]->x, face->v[2]->x,
            s1, s2, s3, normal, t, collision_epsilon);

    if(checked)
        checked = t >=0 && t<=1 && s1 >= -zxEPSILON && s2 >= -zxEPSILON && s3 >= -zxEPSILON;




    if(checked)
    {

        //std::cout<<"test1: "<<t<<" "<<s1<<" "<<s2<<" "<<s3<<std::endl;
        zxContactPoint::Ptr contact = zxContactPoint::create();
        contact->m_normal = normal;
        contact->m_verts.resize(4);
        contact->m_weights.resize(4);
        contact->m_verts[0] = vert;
        contact->m_verts[1] = face->v[0];
        contact->m_verts[2] = face->v[1];
        contact->m_verts[3] = face->v[2];

        contact->m_weights[0] = 1.0;
        contact->m_weights[1] = -s1;
        contact->m_weights[2] = -s2;
        contact->m_weights[3] = -s3;
        contact->m_thickness = m_margin * 0.9;
        m_contacts.push_back(contact);

    }

}

void zxFaceFace_CCD::do_EE_Test(zxCollisionMesh::Edge::Ptr edge0, zxCollisionMesh::Edge::Ptr edge1)
{



    real s0,s1,s2,s3;
    vec3d normal;
    real t;
    real collision_epsilon = zxEPSILON;

    bool checked = zx_check_edge_edge_collision(
                edge0->v[0]->xp, edge0->v[1]->xp, edge1->v[0]->xp, edge1->v[1]->xp,
            edge0->v[0]->x, edge0->v[1]->x, edge1->v[0]->x, edge1->v[1]->x,
            s0, s2, normal, t, collision_epsilon);


    if(checked)
        checked = t >=0 && t<=1 && s0 > zxEPSILON && s0 < 1 - zxEPSILON && s2 > zxEPSILON && s2 < 1 - zxEPSILON;

    s1 = 1 - s0;
    s3 = 1 - s2;

    if(checked)
    {

        zxContactPoint::Ptr contact = zxContactPoint::create();
        contact->m_normal = normal;
        contact->m_verts.resize(4);
        contact->m_weights.resize(4);
        contact->m_verts[0] = edge0->v[0];
        contact->m_verts[1] = edge0->v[1];
        contact->m_verts[2] = edge1->v[0];
        contact->m_verts[3] = edge1->v[1];

        contact->m_weights[0] = s0;
        contact->m_weights[1] = s1;
        contact->m_weights[2] = -s2;
        contact->m_weights[3] = -s3;
        contact->m_thickness = m_margin * 0.9;
        m_contacts.push_back(contact);
    }
}
