#include "zxcglinearsolver.h"

zxCGLinearSolver::zxCGLinearSolver(zxCGMultiplier::Ptr A)
{
    m_A = A;

    m_tol = 1e-6;
    m_maxIter = 100;
}


Eigen::VectorXd zxCGLinearSolver::solver(Eigen::VectorXd& b)
{
    Eigen::VectorXd xk = b;
    xk.setZero();

    Eigen::VectorXd rk,pk;
    Eigen::VectorXd Apk;

    rk = b;
    pk = rk;
    int iter = 0;
    do
    {
        real normRk = rk.squaredNorm();
        Apk = m_A->multiply(pk);
        real pap = pk.dot(Apk);
        real alphaK = normRk / pap;

        if(pap < zxEPSILON)
            break;




        xk += alphaK * pk;
        rk -= alphaK * Apk;

        real normRk1 = rk.squaredNorm();
        if(normRk < m_tol)
            break;

        real beta = normRk1 / normRk;
        pk *= beta;
        pk +=rk;

    }while(iter++ < m_maxIter);

    return xk;

}
