#ifndef ZXBODY_H
#define ZXBODY_H

#include "zxsettings.h"
#include "zxcollisionmesh.h"
#include "zxrendermesh.h"
#include "zxtimestepper.h"

class zxBody
{
    ZX_MAKE_SHARED_MACO_DEFAULT_CREAT(zxBody)
public:
    class Config
    {


    };

public:
    zxBody();


public:
    void                set_render_mesh(zxRenderMesh::Ptr mesh) {m_rendermesh = mesh;}
    void                set_collision_mesh(zxCollisionMesh::Ptr mesh);
    void                set_stepper(zxTimeStepper::Ptr stepper){m_stepper = stepper;}
public:
    zxRenderMesh::Ptr      get_render_mesh(){return m_rendermesh;}
    zxCollisionMesh::Ptr   get_collision_mesh(){return m_colmesh;}
    zxBVHTree::Ptr         get_bvh_tree(){return m_bvh_proxy_tree;}
    zxTimeStepper::Ptr     get_stepper(){return m_stepper;}

    virtual void           do_unconstrained_step(real dt);
    virtual bool           is_static(){return m_stepper == nullptr;}

public:
    virtual void    reset_iac();
    virtual void    compute_iac_response();

protected:
    zxCollisionMesh::Ptr    m_colmesh;
    zxRenderMesh::Ptr       m_rendermesh;
    zxBVHTree::Ptr          m_bvh_proxy_tree, m_bvh_ccd_tree;
    zxTimeStepper::Ptr      m_stepper;

};

#endif // ZXBODY_H
