#ifndef ZXCUBATUREMODELFEM_H
#define ZXCUBATUREMODELFEM_H

#include "zxforcemodel.h"
#include "zxcubaturemodel.h"

class zxCubatureModelFEM : public zxCubatureModel
{
    ZX_MAKE_SHARED_MACO(zxCubatureModelFEM)
public:
    zxCubatureModelFEM(zxForceModel::Ptr forcemodel,Eigen::MatrixXd& U,Eigen::VectorXd& Lambda);

    static Ptr create(zxForceModel::Ptr forcemodel,Eigen::MatrixXd& U,Eigen::VectorXd& Lambda)
    {
        return Ptr (new zxCubatureModelFEM(forcemodel,U,Lambda));
    }

    virtual Eigen::VectorXd     evaluateY(Eigen::VectorXd& X);
    virtual Eigen::VectorXd     get_element_column(int el) ;
    virtual size_t              get_num_elements();

public:
    zxForceModel::Ptr   m_force_model;
    Eigen::MatrixXd     m_mode_U,m_mode_UTrans;

};

#endif // ZXCUBATUREMODELFEM_H
