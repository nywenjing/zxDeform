#ifndef ZXCLOTHFORCEMODEL_H
#define ZXCLOTHFORCEMODEL_H

#include "zxforcemodel.h"
class zxClothForceModel : public zxForceModel
{
public:
    zxClothForceModel(zxTriangularMesh::Ptr mesh);

public:
    virtual void        updatePosition(const Eigen::VectorXd& u) = 0;
    virtual void        computeForce(Eigen::VectorXd& force) = 0;
    virtual void        computeTangent(zxSparseMatrix& tangentK) = 0;
};

#endif // ZXCLOTHFORCEMODEL_H
