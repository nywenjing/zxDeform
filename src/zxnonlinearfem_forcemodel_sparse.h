#ifndef ZXNONLINEARFEM_FORCEMODEL_SPARSE_H
#define ZXNONLINEARFEM_FORCEMODEL_SPARSE_H

#include "zxforcemodel.h"

class zxNonlinearFEM_ForceModel_Sparse : public zxForceModel
{
    ZX_MAKE_SHARED_MACO(zxNonlinearFEM_ForceModel_Sparse)
public:
    zxNonlinearFEM_ForceModel_Sparse(){}
    zxNonlinearFEM_ForceModel_Sparse(zxSolidMesh::Ptr mesh,zxMaterial::Ptr material,real low_clamp,real upper_clamp);

public:
    static Ptr create(zxSolidMesh::Ptr mesh,zxMaterial::Ptr material,real low_clamp = 0.2,real upper_clamp = 1.8)
    {
        return Ptr(new zxNonlinearFEM_ForceModel_Sparse(mesh,material,low_clamp,upper_clamp));
    }

public:
    virtual void        updatePosition(const Eigen::VectorXd& u);
    virtual void        computeForce(Eigen::VectorXd& force);
    virtual void        computeTangent(zxSparseMatrix& tangentK);

    virtual void        set_svd_lowerBound(real lb){ m_low_svd_diag = lb;}
    virtual void        set_svd_upperBound(real ub){ m_upper_svd_diag = ub;}

public:

    real        m_low_svd_diag;
    real        m_upper_svd_diag;


};

#endif // ZXNONLINEARFEM_FORCEMODEL_SPARSE_H
