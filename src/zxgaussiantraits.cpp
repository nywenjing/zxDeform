#include "zxgaussiantraits.h"

template<size_t NN,size_t NI>
void zxGaussianTrait_Interface<NN,NI>::init()
{
    std::vector<real> H(NN);
    for(size_t ni = 0; ni < NI; ni++)
    {
        shape_func(H.data(),m_grst(ni,0),m_grst(ni,1),m_grst(ni,2));

        for(size_t nn = 0; nn < NN; nn++)
            m_shape_func(nn,ni) = H[nn];
    }

    std::vector<real> Hr(NN),Hs(NN),Ht(NN);

    for(size_t ni = 0; ni < NI; ni++)
    {
        shape_deriv(Hr.data(),Hs.data(),Ht.data(),m_grst(ni,0),m_grst(ni,1),m_grst(ni,2));

        for(size_t nn = 0; nn < NN; nn++)
        {
            m_dgdr(nn,ni) = Hr[nn];
            m_dgds(nn,ni) = Hs[nn];
            m_dgdt(nn,ni) = Ht[nn];

        }
    }


}

zxTetrahedron4G1Trait::zxTetrahedron4G1Trait()
{
    m_type = C3D4_1;

    m_grst(0,0) = m_grst(0,1) = m_grst(0,2) = 0.25;
    m_gw[0] = 1.0/6;

    init();
}

void zxTetrahedron4G1Trait::shape_func(real* H,real r,real s,real t)
{
    H[0] = 1 - r - s - t;
    H[1] = r;
    H[2] = s;
    H[3] = t;

}

void zxTetrahedron4G1Trait::shape_deriv(real* Hr,real* Hs,real* Ht,real r, real s, real t)
{
    Hr[0] = -1; Hs[0] = -1; Ht[0] = -1;
    Hr[1] =  1;	Hs[1] =  0; Ht[1] =  0;
    Hr[2] =  0;	Hs[2] =  1; Ht[2] =  0;
    Hr[3] =  0;	Hs[3] =  0; Ht[3] =  1;

}

zxTetrahedron10G4Trait::zxTetrahedron10G4Trait()
{
    m_type = C3D10_4;
    // integration point coordinates
    const double a = 0.58541020;
    const double b = 0.13819660;
    const double w = 0.25 / 6.0;
    m_grst(0,0) = a; m_grst(0,1) = b; m_grst(0,2) = b; m_gw[ 0] = w;
    m_grst(1,0) = b; m_grst(1,1) = a; m_grst(1,2) = b; m_gw[ 0] = w;
    m_grst(2,0) = b; m_grst(2,1) = b; m_grst(2,2) = a; m_gw[ 0] = w;
    m_grst(3,0) = b; m_grst(3,1) = b; m_grst(3,2) = b; m_gw[ 0] = w;

    init();
}

void zxTetrahedron10G4Trait::shape_func(real* H,real r,real s,real t)
{
    double r1 = 1.0 - r - s - t;
    double r2 = r;
    double r3 = s;
    double r4 = t;

    H[0] = r1*(2.0*r1 - 1.0);
    H[1] = r2*(2.0*r2 - 1.0);
    H[2] = r3*(2.0*r3 - 1.0);
    H[3] = r4*(2.0*r4 - 1.0);
    H[4] = 4.0*r1*r2;
    H[5] = 4.0*r2*r3;
    H[6] = 4.0*r3*r1;
    H[7] = 4.0*r1*r4;
    H[8] = 4.0*r2*r4;
    H[9] = 4.0*r3*r4;

}

void zxTetrahedron10G4Trait::shape_deriv(real* Hr,real* Hs,real* Ht,real r, real s, real t)
{
    Hr[0] = -3.0 + 4.0*r + 4.0*(s + t);
    Hr[1] =  4.0*r - 1.0;
    Hr[2] =  0.0;
    Hr[3] =  0.0;
    Hr[4] =  4.0 - 8.0*r - 4.0*(s + t);
    Hr[5] =  4.0*s;
    Hr[6] = -4.0*s;
    Hr[7] = -4.0*t;
    Hr[8] =  4.0*t;
    Hr[9] =  0.0;

    Hs[0] = -3.0 + 4.0*s + 4.0*(r + t);
    Hs[1] =  0.0;
    Hs[2] =  4.0*s - 1.0;
    Hs[3] =  0.0;
    Hs[4] = -4.0*r;
    Hs[5] =  4.0*r;
    Hs[6] =  4.0 - 8.0*s - 4.0*(r + t);
    Hs[7] = -4.0*t;
    Hs[8] =  0.0;
    Hs[9] =  4.0*t;

    Ht[0] = -3.0 + 4.0*t + 4.0*(r + s);
    Ht[1] =  0.0;
    Ht[2] =  0.0;
    Ht[3] =  4.0*t - 1.0;
    Ht[4] = -4.0*r;
    Ht[5] =  0.0;
    Ht[6] = -4.0*s;
    Ht[7] =  4.0 - 8.0*t - 4.0*(r + s);
    Ht[8] =  4.0*r;
    Ht[9] =  4.0*s;

}


zxHexahedronG8Trait::zxHexahedronG8Trait()
{
    m_type = C3D8_8;


    const double a = 1.0 / sqrt(3.0);

    m_grst(0,0) = -a; m_grst(0,1) = -a; m_grst(0,2) = -a; m_gw[0] = 1.0/8;
    m_grst(1,0) =  a; m_grst(1,1) = -a; m_grst(1,2) = -a; m_gw[1] = 1.0/8;
    m_grst(2,0) =  a; m_grst(2,1) =  a; m_grst(2,2) = -a; m_gw[2] = 1.0/8;
    m_grst(3,0) = -a; m_grst(3,1) =  a; m_grst(3,2) = -a; m_gw[3] = 1.0/8;
    m_grst(4,0) = -a; m_grst(4,1) = -a; m_grst(4,2) =  a; m_gw[4] = 1.0/8;
    m_grst(5,0) =  a; m_grst(5,1) = -a; m_grst(5,2) =  a; m_gw[5] = 1.0/8;
    m_grst(6,0) =  a; m_grst(6,1) =  a; m_grst(6,2) =  a; m_gw[6] = 1.0/8;
    m_grst(7,0) = -a; m_grst(7,1) =  a; m_grst(7,2) =  a; m_gw[7] = 1.0/8;

    init();
}

void zxHexahedronG8Trait::shape_func(real* H,real r,real s,real t)
{
    H[0] = 0.125*(1 - r)*(1 - s)*(1 - t);
    H[1] = 0.125*(1 + r)*(1 - s)*(1 - t);
    H[2] = 0.125*(1 + r)*(1 + s)*(1 - t);
    H[3] = 0.125*(1 - r)*(1 + s)*(1 - t);
    H[4] = 0.125*(1 - r)*(1 - s)*(1 + t);
    H[5] = 0.125*(1 + r)*(1 - s)*(1 + t);
    H[6] = 0.125*(1 + r)*(1 + s)*(1 + t);
    H[7] = 0.125*(1 - r)*(1 + s)*(1 + t);

}

void zxHexahedronG8Trait::shape_deriv(real* Hr,real* Hs,real* Ht,real r, real s, real t)
{
    Hr[0] = -0.125*(1 - s)*(1 - t);
    Hr[1] =  0.125*(1 - s)*(1 - t);
    Hr[2] =  0.125*(1 + s)*(1 - t);
    Hr[3] = -0.125*(1 + s)*(1 - t);
    Hr[4] = -0.125*(1 - s)*(1 + t);
    Hr[5] =  0.125*(1 - s)*(1 + t);
    Hr[6] =  0.125*(1 + s)*(1 + t);
    Hr[7] = -0.125*(1 + s)*(1 + t);

    Hs[0] = -0.125*(1 - r)*(1 - t);
    Hs[1] = -0.125*(1 + r)*(1 - t);
    Hs[2] =  0.125*(1 + r)*(1 - t);
    Hs[3] =  0.125*(1 - r)*(1 - t);
    Hs[4] = -0.125*(1 - r)*(1 + t);
    Hs[5] = -0.125*(1 + r)*(1 + t);
    Hs[6] =  0.125*(1 + r)*(1 + t);
    Hs[7] =  0.125*(1 - r)*(1 + t);

    Ht[0] = -0.125*(1 - r)*(1 - s);
    Ht[1] = -0.125*(1 + r)*(1 - s);
    Ht[2] = -0.125*(1 + r)*(1 + s);
    Ht[3] = -0.125*(1 - r)*(1 + s);
    Ht[4] =  0.125*(1 - r)*(1 - s);
    Ht[5] =  0.125*(1 + r)*(1 - s);
    Ht[6] =  0.125*(1 + r)*(1 + s);
    Ht[7] =  0.125*(1 - r)*(1 + s);

}

zxTriangle3G3Trait::zxTriangle3G3Trait()
{
    m_type = S3_3;

    real a = 1.0/6.0;
    real b = 2.0/3.0;

    m_grst(0,0) = a; m_grst(0,1) = a; m_grst(0,2) = 0.0;
    m_grst(1,0) = b; m_grst(1,1) = a; m_grst(1,2) = 0.0;
    m_grst(2,0) = a; m_grst(2,1) = b; m_grst(2,2) = 0.0;

    m_gw[0] = m_gw[1] = m_gw[2] = 1.0/6;

    init();
}

void zxTriangle3G3Trait::shape_func(real* H,real r,real s,real t)
{
    H[0] = 1.0 - r - s;
    H[1] = r;
    H[2] = s;

}

void zxTriangle3G3Trait::shape_deriv(real* Hr,real* Hs,real* Ht,real r, real s, real t)
{
    Hr[0] = -1; Hs[0] = -1;
    Hr[1] =  1; Hs[1] =  0;
    Hr[2] =  0; Hs[2] =  1;

}

zxTriangle6G3Trait::zxTriangle6G3Trait()
{
    m_type = S6_3;

    real a = 1.0/6.0;
    real b = 2.0/3.0;

    m_grst(0,0) = a; m_grst(0,1) = a; m_grst(0,2) = 0.0;
    m_grst(1,0) = b; m_grst(1,1) = a; m_grst(1,2) = 0.0;
    m_grst(2,0) = a; m_grst(2,1) = b; m_grst(2,2) = 0.0;

    m_gw[0] = m_gw[1] = m_gw[2] = 1.0/6;

    init();
}

void zxTriangle6G3Trait::shape_func(real* H,real r,real s,real t)
{
    double r1 = 1.0 - r - s;
    double r2 = r;
    double r3 = s;

    H[0] = r1*(2.0*r1 - 1.0);
    H[1] = r2*(2.0*r2 - 1.0);
    H[2] = r3*(2.0*r3 - 1.0);
    H[3] = 4.0*r1*r2;
    H[4] = 4.0*r2*r3;
    H[5] = 4.0*r3*r1;

}

void zxTriangle6G3Trait::shape_deriv(real* Hr,real* Hs,real* Ht,real r, real s, real t)
{
    Hr[0] = -3.0 + 4.0*r + 4.0*s;
    Hr[1] =  4.0*r - 1.0;
    Hr[2] =  0.0;
    Hr[3] =  4.0 - 8.0*r - 4.0*s;
    Hr[4] =  4.0*s;
    Hr[5] = -4.0*s;

    Hs[0] = -3.0 + 4.0*s + 4.0*r;
    Hs[1] =  0.0;
    Hs[2] =  4.0*s - 1.0;
    Hs[3] = -4.0*r;
    Hs[4] =  4.0*r;
    Hs[5] =  4.0 - 8.0*s - 4.0*r;

}

zxGaussianFactory::zxGaussianFactory()
{
    m_traits_group[zxTetrahedron4G1Trait::Singleton().m_type] = &(zxTetrahedron4G1Trait::Singleton());
    m_traits_group[zxTetrahedron10G4Trait::Singleton().m_type] = &(zxTetrahedron10G4Trait::Singleton());
    m_traits_group[zxHexahedronG8Trait::Singleton().m_type] = &(zxHexahedronG8Trait::Singleton());
    m_traits_group[zxTriangle3G3Trait::Singleton().m_type] = &(zxTriangle3G3Trait::Singleton());
    m_traits_group[zxTriangle6G3Trait::Singleton().m_type] = &(zxTriangle6G3Trait::Singleton());
}

template<size_t NN,size_t NI>
zxGaussianTrait_Interface<NN,NI>::zxGaussianTrait_Interface()
{
    m_gw = Eigen::VectorXd(NI);
    m_shape_func.resize(NN,NI);
    m_dgdr.resize(NN,NI);
    m_dgds.resize(NN,NI);
    m_dgdt.resize(NN,NI);
    m_grst.resize(NI,3);
}
