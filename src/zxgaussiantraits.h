#ifndef ZXGAUSSIANTRAITS_H
#define ZXGAUSSIANTRAITS_H

#include "zxsettings.h"

class zxGaussianTrait
{
public:
    typedef std::shared_ptr<zxGaussianTrait> Ptr;
    typedef std::shared_ptr<zxGaussianTrait const> ConstPtr;
public:
    enum Type
    {
        INVALID = -1,
        C3D4_1 = 0,
        C3D10_4,
        C3D8_8,
        NumberType
    };

public:
    zxGaussianTrait(){ m_type = INVALID;}

    virtual size_t get_num_gaussian_points() = 0;
    virtual real&   get_shape_fun(int n_id,int g_id) = 0;
    virtual real&   get_first_derive_r(int n_id,int g_id) = 0;
    virtual real&   get_first_derive_s(int n_id,int g_id) = 0;
    virtual real&   get_first_derive_t(int n_id,int g_id) = 0;
    virtual real&   get_gaussian_weight(int g_id) = 0;

public:
    Type    m_type;

};

template<size_t NN,size_t NI>
class zxGaussianTrait_Interface : public zxGaussianTrait
{
public:
    zxGaussianTrait_Interface();
    virtual size_t get_num_gaussian_points() {return NI;}
    virtual real&   get_shape_fun(int n_id,int g_id) {return m_shape_func(n_id,g_id);}
    virtual real&   get_first_derive_r(int n_id,int g_id) {return m_dgdr(n_id,g_id);}
    virtual real&   get_first_derive_s(int n_id,int g_id) {return m_dgds(n_id,g_id);}
    virtual real&   get_first_derive_t(int n_id,int g_id) {return m_dgdt(n_id,g_id);}
    virtual real&   get_gaussian_weight(int g_id) {return m_gw[g_id];}

protected:
    virtual void shape_func(real* H,real r,real s,real t) = 0;
    virtual void shape_deriv(real* Hr,real* Hs,real* Ht,real r, real s, real t) = 0;
    virtual void init();

public:
    Eigen::VectorXd m_gw;
    Eigen::MatrixXd m_shape_func;
    Eigen::MatrixXd m_dgdr;
    Eigen::MatrixXd m_dgds;
    Eigen::MatrixXd m_dgdt;
    Eigen::MatrixXd m_grst;
};

template<size_t NI>
class zxTetrahedron4Trait : public zxGaussianTrait_Interface<4,NI>
{

};

template<size_t NI>
class zxTetrahedron10Trait : public zxGaussianTrait_Interface<10,NI>
{

};

class zxTetrahedron4G1Trait : public zxTetrahedron4Trait<1>
{
private:
    zxTetrahedron4G1Trait();
public:
    typedef std::shared_ptr<zxTetrahedron4G1Trait> Ptr;
    typedef std::shared_ptr<zxTetrahedron4G1Trait const> ConstPtr;
public:
    static zxTetrahedron4G1Trait& Singleton(){static zxTetrahedron4G1Trait instance; return instance;}
public:
    virtual void shape_func(real* H,real r,real s,real t);
    virtual void shape_deriv(real* Hr,real* Hs,real* Ht,real r, real s, real t);

};

class zxTetrahedron10G4Trait : public zxTetrahedron10Trait<4>
{

};

class zxGaussianFactory
{
private:
    zxGaussianFactory();
public:
    static zxGaussianFactory& Singleton(){static zxGaussianFactory instance; return instance;}
public:
    size_t get_num_gaussian_points(zxGaussianTrait::Type type){  return Singleton().get_traint(type)->get_num_gaussian_points(); }
    real&   get_shape_fun(zxGaussianTrait::Type type, int n_id,int g_id) { return Singleton().get_traint(type)->get_shape_fun(n_id,g_id);}
    real&   get_first_derive_r(zxGaussianTrait::Type type, int n_id,int g_id) {return Singleton().get_traint(type)->get_first_derive_r(n_id,g_id);}
    real&   get_first_derive_s(zxGaussianTrait::Type type, int n_id,int g_id) {return Singleton().get_traint(type)->get_first_derive_s(n_id,g_id);}
    real&   get_first_derive_t(zxGaussianTrait::Type type, int n_id,int g_id) {return Singleton().get_traint(type)->get_first_derive_t(n_id,g_id);}
    real&   get_gaussian_weight(zxGaussianTrait::Type type, int g_id) {return Singleton().get_traint(type)->get_gaussian_weight(g_id);}
protected:
    zxGaussianTrait* get_traint(zxGaussianTrait::Type type){return m_traits_group[type];}

public:
    zxGaussianTrait* m_traits_group[zxGaussianTrait::NumberType];
};

#endif // ZXGAUSSIANTRAITS_H
