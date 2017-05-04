#ifndef ZXCOLLISIONALGORITHM_H
#define ZXCOLLISIONALGORITHM_H

#include "zxbvh.h"
#include "zxcollisionmesh.h"
#include "zxcontactpoint.h"
class zxCollisionAlgorithm
{
public:
    zxCollisionAlgorithm();
};

class zxFaceFace_Collider : public zxBVHCollider
{
    ZX_MAKE_SHARED_MACO(zxFaceFace_Collider)
public:
    zxFaceFace_Collider(){ }
public:
    void process_collision(zxBVHNode* node0, zxBVHNode* node1);

public:
    virtual void do_VF_Test(zxCollisionMesh::Vert::Ptr vert, zxCollisionMesh::Face::Ptr face) = 0;
    virtual void do_EE_Test(zxCollisionMesh::Edge::Ptr e0, zxCollisionMesh::Edge::Ptr e1) = 0;

    std::list<zxContactPoint::Ptr>&  get_contacts(){return m_contacts;}

protected:
    std::list<zxContactPoint::Ptr>  m_contacts;


};

class zxFaceFace_Proximity : public zxFaceFace_Collider
{
    ZX_MAKE_SHARED_MACO_DEFAULT_CREAT(zxFaceFace_Proximity)
public:
    virtual void do_VF_Test(zxCollisionMesh::Vert::Ptr vert, zxCollisionMesh::Face::Ptr face);
    virtual void do_EE_Test(zxCollisionMesh::Edge::Ptr e0, zxCollisionMesh::Edge::Ptr e1);

};

class zxFaceFace_CCD : public zxFaceFace_Collider
{
    ZX_MAKE_SHARED_MACO_DEFAULT_CREAT(zxFaceFace_CCD)
public:
    virtual void do_VF_Test(zxCollisionMesh::Vert::Ptr vert, zxCollisionMesh::Face::Ptr face);
    virtual void do_EE_Test(zxCollisionMesh::Edge::Ptr e0, zxCollisionMesh::Edge::Ptr e1);
};

#endif // ZXCOLLISIONALGORITHM_H
