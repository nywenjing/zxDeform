#ifndef ZXMATERIAL_H
#define ZXMATERIAL_H

#include "zxsettings.h"

class zxMaterial
{
    ZX_MAKE_SHARED_MACO(zxMaterial)
public:
    zxMaterial();
    zxMaterial(real density){ m_density = density;}
public:
    real get_density(){return m_density;}

    virtual void compute_gradient(const vec3d& strainInvariants,vec3d& gradient) const = 0;
    virtual void compute_hessian(const vec3d& strainInvariants, real hessian[6])const = 0;
protected:
    real    m_density;
};

#endif // ZXMATERIAL_H
