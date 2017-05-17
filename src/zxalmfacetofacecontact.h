#ifndef ZXALMFACETOFACECONTACT_H
#define ZXALMFACETOFACECONTACT_H

#include "zxmesh.h"
#include "zxcollisionmesh.h"
#include "zxalmcontactinterface.h"
#include "zxgaussiantraits.h"
class zxALMFaceToFaceContact : public zxALMContactInterface
{
    ZX_MAKE_SHARED_MACO(zxALMFaceToFaceContact)
public:
    zxALMFaceToFaceContact(zxCollisionMesh::Ptr mastersurf,zxCollisionMesh::Ptr slavesurf,real eps = 1e7);

    static Ptr create(zxCollisionMesh::Ptr mastersurf,zxCollisionMesh::Ptr slavesurf,real eps = 1e7){return Ptr (new zxALMFaceToFaceContact(mastersurf,slavesurf,eps));}

public:

    virtual void                                update_contact();
    virtual void                                compute_contact_force(Eigen::VectorXd& R);
    virtual std::list<Eigen::Triplet<real>>     get_stiffness_triplets();
    virtual bool                                alm_augment();
    virtual void                                compute_stiffness(Eigen::SparseMatrix<real>& tangentStiff);

public:
    zxCollisionMesh::Ptr    m_masterSurface,m_slaveSurface;
    zxGaussianTrait::Type   m_gaussian_type;

};

#endif // ZXALMFACETOFACECONTACT_H
