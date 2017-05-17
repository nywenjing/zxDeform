#ifndef ZXALMCONTACTINTERFACE_H
#define ZXALMCONTACTINTERFACE_H

#include "zxsettings.h"
#include "zxcollisionmesh.h"
class zxALMContactInterface
{
    ZX_MAKE_SHARED_MACO(zxALMContactInterface)
public:
    class ALMData
    {
    public:
        ALMData()
        {
            m_lambda = 0;
            m_gap = 0;
            m_normal = vec3d::Zero();
            m_master_element = nullptr;
        }

    public:
        real                        m_lambda;
        real                        m_gap;
        vec3d                       m_normal;
        zxCollisionMesh::Face*      m_master_element;
        size_t                      m_id;
        vec3d                       m_x;
        real                        m_r,m_s;
    };

public:
    zxALMContactInterface(real eps);
public:
    virtual void                                update_contact() = 0;
    virtual void                                compute_contact_force(Eigen::VectorXd& R) = 0;
    virtual std::list<Eigen::Triplet<real>>   get_stiffness_triplets() = 0;
    virtual void                                compute_stiffness(Eigen::SparseMatrix<real>& tangentStiff) = 0;
    virtual bool                                alm_augment() = 0;


protected:
    std::vector<std::vector<ALMData>>   m_almData;
    size_t                              m_num_gaussian;
    real                                m_eps;
    real                                m_alm_tol;
};

#endif // ZXALMCONTACTINTERFACE_H
