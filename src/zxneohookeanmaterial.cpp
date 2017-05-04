#include "zxneohookeanmaterial.h"

zxNeoHookeanMaterial::zxNeoHookeanMaterial(real Y_m,real possion_ratio,real density)
    :zxMaterial(density)
{
    mLambda = possion_ratio * Y_m / ((1.0 + possion_ratio) * (1.0 - 2.0 * possion_ratio));
    mMu = Y_m / (2.0 * (1.0 + possion_ratio));
}

void zxNeoHookeanMaterial::compute_gradient(const vec3d& strainInvariants,vec3d& gradient) const
{
    double IIIC = strainInvariants[2];
    gradient[0] = 0.5 * mMu;
    gradient[1] = 0.0;
    gradient[2] = (-0.5 * mMu + 0.25 * mLambda * log(IIIC)) / IIIC;
}

void zxNeoHookeanMaterial::compute_hessian(const vec3d& strainInvariants, real hessian[6]) const
{
    for(int i = 0; i < 6; i++)
        hessian[i] = 0.0;
    double IIIC = strainInvariants[2];
    hessian[5] = (0.25 * mLambda + 0.5 * mMu - 0.25 * mLambda * log(IIIC)) / (IIIC * IIIC);
}
