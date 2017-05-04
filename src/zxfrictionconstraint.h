#ifndef ZXFRICTIONCONSTRAINT_H
#define ZXFRICTIONCONSTRAINT_H

#include "zxcontactconstraint.h"
class zxFrictionConstraint : public zxLCPConstraint
{
    ZX_MAKE_SHARED_MACO(zxFrictionConstraint)
public:
    zxFrictionConstraint();
    zxFrictionConstraint(zxContactConstraint::Ptr cc,vec3d dir);

    static Ptr create(zxContactConstraint::Ptr cc,vec3d dir){return Ptr (new zxFrictionConstraint(cc,dir));}

public:
    virtual void do_projected_gs_solve();
    virtual void set_friction(real friction){m_friction = friction;}

protected:
    zxContactConstraint::Ptr    m_contact;
    real                        m_friction;
};

#endif // ZXFRICTIONCONSTRAINT_H
