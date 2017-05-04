#include "zxbody.h"

zxBody::zxBody()
{
    m_colmesh = nullptr;
    m_rendermesh = nullptr;
    m_bvh_proxy_tree = nullptr;
    m_bvh_ccd_tree = nullptr;
    m_stepper = nullptr;

}

void zxBody::do_unconstrained_step(real dt)
{
    m_stepper->do_step(dt);
    m_stepper->update_mesh();
}

void zxBody::reset_iac()
{
    m_stepper->reset_iac();
}

void zxBody::compute_iac_response()
{
    m_stepper->compute_iac_response();
}

void zxBody::set_collision_mesh(zxCollisionMesh::Ptr mesh)
{
    m_colmesh = mesh;
    std::vector<zxCollisionMesh::Face::Ptr>& faces = m_colmesh->get_faces();
    std::vector<zxAABBData::Ptr> data;
    data.insert(data.end(),faces.begin(),faces.end());
    m_bvh_proxy_tree = zxBVHTree::create(data);
}
