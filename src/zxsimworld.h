#ifndef ZXSIMWORLD_H
#define ZXSIMWORLD_H

#include "zxbody.h"
#include "zxlcpconstraint.h"
#include "zxcontactpoint.h"

class zxSimWorld
{
    ZX_MAKE_SHARED_MACO_DEFAULT_CREAT(zxSimWorld)
public:
    zxSimWorld();

public:
    size_t      get_num_bodies(){return m_bodies.size();}
    zxBody::Ptr get_body(size_t id){return m_bodies.at(id);}

public:
    void        add_body(zxBody::Ptr body){m_bodies.push_back(body);}

public:
    void        do_time_step();
public:
    void        set_margin(real margin){ m_margin = margin;}
    void        set_friction(real friction){ m_friction = friction;}
protected:
    void        do_unconstrained_step();
    size_t      do_ccd();
    size_t      do_proxy();
    void        do_iac();
protected:
    std::vector<zxBody::Ptr>            m_bodies;
    std::list<zxLCPConstraint::Ptr>     m_dynamic_lcps;
    std::vector<zxLCPConstraint::Ptr>   m_static_lcps;
    real                                m_dt;
    real                                m_margin;
    real                                m_friction;

public:
    std::vector<zxContactPoint::Ptr>    m_debug_ccd_contact;

};

#endif // ZXSIMWORLD_H
