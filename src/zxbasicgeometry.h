#ifndef ZXBASICGEOMETRY_H
#define ZXBASICGEOMETRY_H

#include "zxsettings.h"
#include "zxgaussiantraits.h"
#include "zxbvh.h"

class zxMaterialPoint
{
    ZX_MAKE_SHARED_MACO(zxMaterialPoint)
    static Ptr create(){return Ptr (new zxMaterialPoint());}

public:
    zxMaterialPoint()
    {
        m_dPsdF_diag.resize(9,9);
        m_dFdr.resize(9,12);
        m_dFdr_trans.resize(12,9);
    }

public:
    mat3d m_defgrad;
    mat3d m_invJac0;
    mat3d m_Pstress;
    mat3d m_svd_U;
    mat3d m_svd_V;
    vec3d m_svd_diag;
    vec3d m_diag_P;
    real  m_detJac0;
    Eigen::MatrixXd m_dPsdF_diag;
    Eigen::MatrixXd m_dFdr,m_dFdr_trans;

};

enum zxBC
{
    zxFree = 0,
    zxFixed,
    zxPrescribed
};

class zxNode
{
ZX_MAKE_SHARED_MACO(zxNode)
public:
    zxNode()
    {
        rp.setZero();
        rt.setZero();
        r0.setZero();

        m_bc[0] = m_bc[1] = m_bc[2] = zxFree;
        m_dof_id[0] = m_dof_id[1] = m_dof_id[2] = -1;

        reset_iac();
    }

static Ptr create(){return Ptr (new zxNode());}

public:
    vec3d rp,rt,r0;
    int   m_id;

public:
    vec3d rl;

public:
    void reset_iac(){ m_isLcp = false; m_lcp_res_imp = m_lcp_dx = m_lcp_pcg_dx = m_lcp_imp = vec3d::Zero(); m_lcp_invA.setZero();}
    bool    m_isLcp;
    vec3d   m_lcp_dx,m_lcp_pcg_dx,m_lcp_imp,m_lcp_res_imp;
    mat3d   m_lcp_invA;

public:

    zxBC    m_bc[3];
    int     m_dof_id[3];

};

class zxElement
{
ZX_MAKE_SHARED_MACO(zxElement)
public:
    size_t         get_num_nodes(){return m_nodes.size();}
    zxNode::Ptr    get_node(size_t i) {return m_nodes.at(i);}

public:
    void            set_num_nodes(size_t n){m_nodes.resize(n);}
    void            set_node(size_t id, zxNode::Ptr node){m_nodes.at(id) = node;}

public:
    virtual real    get_volume(){return 0;}

public:
    virtual size_t get_num_gaussian_points() {return zxGaussianFactory::Singleton().get_num_gaussian_points(m_type);}
    virtual real&   get_shape_fun(int n_id,int g_id) {return zxGaussianFactory::Singleton().get_shape_fun(m_type,n_id,g_id);}
    virtual real&   get_first_derive_r(int n_id,int g_id) {return zxGaussianFactory::Singleton().get_first_derive_r(m_type,n_id,g_id);}
    virtual real&   get_first_derive_s(int n_id,int g_id) {return zxGaussianFactory::Singleton().get_first_derive_s(m_type,n_id,g_id);}
    virtual real&   get_first_derive_t(int n_id,int g_id) {return zxGaussianFactory::Singleton().get_first_derive_t(m_type,n_id,g_id);}
    virtual real&   get_gaussian_weight(int g_id) {return zxGaussianFactory::Singleton().get_gaussian_weight(m_type,g_id);}
    zxMaterialPoint::Ptr get_material_point(size_t m_id) {return m_material_points.at(m_id);}
public:
    std::vector<zxNode::Ptr> m_nodes;
    std::vector<zxMaterialPoint::Ptr> m_material_points;
    zxGaussianTrait::Type             m_type;
    size_t                            m_id;
};

class zxSurfaceElement : public zxElement
{
public:
    virtual real get_area() = 0;

};

class zxSolidElement : public zxElement
{
ZX_MAKE_SHARED_MACO(zxElement)
public:
    virtual real get_volume() = 0;
    virtual void init_material_points();

};

class zxTriangle : public zxSurfaceElement
{
ZX_MAKE_SHARED_MACO(zxTriangle)

public:
    real get_area();



};

class zxTetrahedron : public zxSolidElement
{
ZX_MAKE_SHARED_MACO(zxTetrahedron)

static Ptr create(){return Ptr (new zxTetrahedron());}
public:
    real get_volume();

};

#endif // ZXBASICGEOMETRY_H
