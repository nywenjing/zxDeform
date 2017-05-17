#ifndef ZXALMSIMULATOR_H
#define ZXALMSIMULATOR_H

#include "zxforcemodel.h"
#include "zxalmcontactinterface.h"

//#include "Eigen/PardisoSupport"

class zxNodalForce
{
    ZX_MAKE_SHARED_MACO_DEFAULT_CREAT(zxNodalForce)
public:
    zxNode* m_node;
    vec3d   m_force;
};

class zxALMSimulator
{
    ZX_MAKE_SHARED_MACO_DEFAULT_CREAT(zxALMSimulator)
public:
    zxALMSimulator();

public:
    void    add_forcemodel(zxForceModel::Ptr forcemodel);
    void    add_contactInterface(zxALMContactInterface::Ptr contactInterface);
    void    add_nodal_force(zxNodalForce::Ptr nforce){ m_nodal_forces.push_back(nforce);}

public:
    void    init_simulator();
    void    do_simulate();
protected:
    void    init_step();
    void    update_deformation(Eigen::VectorXd& u);
    void    update_contact();
    void    update_stress();
    void    compute_residual(Eigen::VectorXd& R);
    void    compute_stiffness_matrix();
    void    reform_stiffness_topology();

    real    do_line_search(real s);
    bool    alm_augment();



protected:
    std::vector<zxNode::Ptr>        m_nodes;
    std::vector<zxForceModel::Ptr>  m_forceModels;
    std::vector<zxALMContactInterface::Ptr> m_contactInterfaces;
    size_t                          m_numDofs;

protected:
    std::vector<zxNodalForce::Ptr>       m_nodal_forces;

protected:
    Eigen::VectorXd                                     m_ui,m_Ui,m_Ut,m_R0,m_R1,m_Fd;
    zxSparseMatrix                                      m_stiffness_matrix;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<real>>        m_linearSolver;
    std::vector<bool>                                    m_is_prescribed_dofs;

protected:
    real    m_residual_tol;
    real    m_disp_tol;
    real    m_energy_tol;
};

#endif // ZXALMSIMULATOR_H
