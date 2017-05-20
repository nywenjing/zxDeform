#ifndef ZXTIMESTEPPER_H
#define ZXTIMESTEPPER_H

#include "zxsettings.h"
#include "zxsparsematrix.h"
#include "zxforcemodel.h"

class zxTimeStepper
{
public:
    typedef std::shared_ptr<zxTimeStepper> Ptr;
    typedef std::shared_ptr<zxTimeStepper const> ConstPtr;
public:
    zxTimeStepper();
    zxTimeStepper(zxForceModel::Ptr forcemodel);

    virtual ~zxTimeStepper(){}

    //static Ptr create(zxForceModel::Ptr forcemodel){ return Ptr (new zxTimeStepper(forcemodel));}



public:
    virtual void do_step(real dt) = 0;
    virtual void compute_iac_response() = 0;
    virtual void reset_iac() = 0;
    virtual void update_mesh();
    virtual size_t  get_dim() = 0;
    virtual zxForceModel::Ptr   get_force_model(){return m_forceModel;}

public:
    virtual void do_cg_alm_multi(double* Ab,double* b) = 0;
    virtual void do_project_alm_ab(double *Ab) = 0;
    virtual void do_project_alm_impulse(double* b) = 0;
    virtual void post_update_alm_constraint(double* x) = 0;


protected:
    Eigen::VectorXd m_pos,m_pos0,m_vel, m_vel0,m_dpos, m_ext_force,m_int_force,m_tem_pos;

protected:
    zxForceModel::Ptr   m_forceModel;
    size_t              m_dim;
    size_t              m_unconstrained_dim;
    std::vector<size_t> m_fixedDof;

protected:
    zxSparseMatrix  m_mass_spmat;
    zxSparseMatrix  m_stiff_spmat;
    zxSparseMatrix  m_sys_spmat;
    Eigen::VectorXd m_residual;

protected:
    Eigen::MatrixXd m_mass_dnmat,m_stfiff_dnmat,m_sys_dnmat;

    real            m_dt;
};

#endif // ZXTIMESTEPPER_H
