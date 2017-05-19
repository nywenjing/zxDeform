#include "zxcubaturemodel.h"
#include "mersennetwister.h"

zxCubatureModel::zxCubatureModel()
{

}


void zxCubatureModel::generateSample(int numSample, real magnitude)
{
    mSampleX.resize(get_X_dim(),numSample);
    mSampleY.resize(get_Y_dim(),numSample);

    int dimX = get_X_dim();


    for (int x = 0; x < dimX; x++)
        if (mLamda[x] < 1e-3)
        {
            mLamda[x] = 0;
        }

    int firstNonZero = 0;
    for (int x = 0; x < dimX; x++)
        if (mLamda[x] > 0.0)
        {
            firstNonZero = x;
            break;
        }
    for (int x = 0; x < firstNonZero; x++)
        mLamda[x] = mLamda[firstNonZero];

    // the user wants the first mode to have the given magnitude
    double _magnitude = magnitude * sqrt(mLamda[0]);

    MERSENNETWISTER _trainingTwister;
    for(int i = 0; i < numSample; i++)
    {
        Eigen::VectorXd X(dimX);
        for(int j = 0; j < dimX; j++)
            X[j] = _magnitude / sqrt(mLamda[j]) / 4.0 * _trainingTwister.randNorm();

        Eigen::VectorXd Y = evaluateY(X);

        mSampleX.col(i) = X;
        mSampleY.col(i) = Y;
    }

    mSampleMagnitudeScale.resize(numSample);

    for(int isample = 0; isample < numSample; isample++)
    {
        double mag = mSampleY.col(isample).norm();
        mag = 1.0/mag;
        mSampleMagnitudeScale[isample] = mag;
        mSampleY.col(isample) *= mag;
    }

}
