#ifndef ZXTIMESTEPPERBACKWARDEULER_H
#define ZXTIMESTEPPERBACKWARDEULER_H

#include "zxtimestepper.h"
#include "zxforcemodel.h"

class zxTimeStepperBackwardEuler : public zxTimeStepper
{
    ZX_MAKE_SHARED_MACO(zxTimeStepperBackwardEuler)

public:
    zxTimeStepperBackwardEuler();
    zxTimeStepperBackwardEuler(zxForceModel::Ptr forcemodel);
    static Ptr create(zxForceModel::Ptr forcemodel){return Ptr (new zxTimeStepperBackwardEuler(forcemodel));}
    public:
        void do_step(real dt);
    protected:
        void do_step_sparse(real dt);
    protected:
        void compute_hessian(zxSparseMatrix& H);
        void compute_residual(Eigen::VectorXd& dres);
        void compute_delta_position(Eigen::VectorXd& dp);
        void update_position(Eigen::VectorXd& q);

public:
        virtual void compute_iac_response();
        virtual void compute_iac_response_gs(Eigen::VectorXd& x,Eigen::VectorXd& b,Eigen::VectorXd &r);
        virtual void reset_iac();

};

#endif // ZXTIMESTEPPERBACKWARDEULER_H
