#ifndef ZXCUBATUREGENERATORNNHTP_H
#define ZXCUBATUREGENERATORNNHTP_H

#include "zxcubaturemodel.h"
#include "mersennetwister.h"
class zxCubatureGeneratorNNHTP
{
    ZX_MAKE_SHARED_MACO(zxCubatureGeneratorNNHTP)
public:
    zxCubatureGeneratorNNHTP(zxCubatureModel::Ptr cmodel);

    static Ptr create(zxCubatureModel::Ptr cmodel){ return Ptr (new zxCubatureGeneratorNNHTP(cmodel));}

public:
    void generateCubatures(int numC);
    void generateSamples(int numSample,double mag);
protected:
    std::vector<int> random_pick_candidates(std::vector<int>& excludes,int candidatePerTry,
                                            int totalCandidates,bool addExcludes);
    std::vector<int>     project(Eigen::VectorXd& w, int toKeep);
    Eigen::VectorXd      compute_element_Quantity(int id);

protected:
    zxCubatureModel::Ptr    m_cubature_model;
    MERSENNETWISTER         m_twister;

protected:
    std::vector<int>    mKeyElements;
    std::vector<double>   mKeyWeights;
};

#endif // ZXCUBATUREGENERATORNNHTP_H
