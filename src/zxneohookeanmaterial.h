#ifndef ZXNEOHOOKEANMATERIAL_H
#define ZXNEOHOOKEANMATERIAL_H

#include "zxmaterial.h"

class zxNeoHookeanMaterial : public zxMaterial
{
    ZX_MAKE_SHARED_MACO(zxNeoHookeanMaterial)
public:
    zxNeoHookeanMaterial(real Y_m,real possion_ratio,real density);

    static Ptr create(real Y_m = 1e6,real possion_ratio = 0.33,real density = 1e3){ return Ptr (new zxNeoHookeanMaterial(Y_m,possion_ratio,density));}

public:
    virtual void compute_gradient(const vec3d& strainInvariants,vec3d& gradient) const;
    virtual void compute_hessian(const vec3d& strainInvariants, real hessian[6]) const;

protected:
    real    mLambda;
    real    mMu;
};

#endif // ZXNEOHOOKEANMATERIAL_H
