#ifndef ZXCONTACTCONSTRAINT_H
#define ZXCONTACTCONSTRAINT_H


#include "zxlcpconstraint.h"
#include "zxcontactpoint.h"

class zxTestConstraint : public zxLCPConstraint
{
public:
    zxTestConstraint(zxContactPoint::Ptr cp)
    {
        cp->m_normal *= 1.0;
    }
};

class zxContactConstraint : public zxLCPConstraint
{
    ZX_MAKE_SHARED_MACO(zxContactConstraint)
public:
    zxContactConstraint(zxContactPoint* cp);
    void setcontact(zxContactPoint::Ptr cp);

    static Ptr create(zxContactPoint* cp){return Ptr(new zxContactConstraint(cp));}
};

#endif // ZXCONTACTCONSTRAINT_H
