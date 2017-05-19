#ifndef ZXCUBATUREMODEL_H
#define ZXCUBATUREMODEL_H

#include "zxsettings.h"
class zxCubatureModel
{
    ZX_MAKE_SHARED_MACO(zxCubatureModel)
public:
    zxCubatureModel();

public:
    size_t                      get_num_samples(){return mSampleX.cols();}
    size_t                      get_X_dim(){return m_dimX;}
    size_t                      get_Y_dim(){return m_dimY;}
    virtual void                generateSample(int numSample,real magnitude);
    Eigen::MatrixXd&            get_sample_X(){return mSampleX;}
    Eigen::MatrixXd&            get_sample_Y(){return mSampleY;}

public:
    virtual Eigen::VectorXd     evaluateY(Eigen::VectorXd& X) = 0;
    virtual Eigen::VectorXd     get_element_column(int el) = 0;
    virtual size_t              get_num_elements() = 0;

protected:
    Eigen::MatrixXd     mSampleX,mSampleY;
    Eigen::VectorXd     mLamda;
    size_t              m_dimX,m_dimY;
    Eigen::VectorXd     mSampleMagnitudeScale;
};

#endif // ZXCUBATUREMODEL_H
