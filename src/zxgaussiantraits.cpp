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

zxGaussianFactory::zxGaussianFactory()
{
    m_traits_group[zxTetrahedron4G1Trait::Singleton().m_type] = &(zxTetrahedron4G1Trait::Singleton());
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
