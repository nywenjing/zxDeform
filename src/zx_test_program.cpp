#include "zxbasicgeometry.h"
#include "zxmesh.h"
#include "zxnonlinearfem_forcemodel_sparse.h"
#include "zxmaterial.h"
#include "zxneohookeanmaterial.h"
#include "zxbody.h"
#include "zxtimestepperbackwardeuler.h"

void test_basic_geometry()
{

    std::string filename = "../zxDeform/data/box.1";
    zxSolidMesh::Ptr t_mesh = zxTetrahedralMesh::create(filename);
    zxMaterial::Ptr material = zxNeoHookeanMaterial::create();
    zxNonlinearFEM_ForceModel_Sparse::Ptr forceModel = zxNonlinearFEM_ForceModel_Sparse::create(t_mesh,material);

    forceModel->getForceDimension();


    zxTimeStepperBackwardEuler::Ptr stepper = zxTimeStepperBackwardEuler::create(forceModel);
    stepper->do_step(0.01);


}
