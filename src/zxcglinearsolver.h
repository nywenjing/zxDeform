#ifndef ZXCGLINEARSOLVER_H
#define ZXCGLINEARSOLVER_H

#include "zxsettings.h"
class zxCGMultiplier
{
    ZX_MAKE_SHARED_MACO(zxCGMultiplier)
public:
    virtual Eigen::VectorXd multiply(Eigen::VectorXd& b) = 0;
};

class zxCGLinearSolver
{
public:
    zxCGLinearSolver(zxCGMultiplier::Ptr A);

public:
    Eigen::VectorXd solver(Eigen::VectorXd& b);

    virtual void set_tolerance(real tol){m_tol = tol;}
    virtual void set_maxIter(int iter){ m_maxIter = iter;}

protected:
    zxCGMultiplier::Ptr m_A;
    real                m_tol;
    int                 m_maxIter;
};

#endif // ZXCGLINEARSOLVER_H
