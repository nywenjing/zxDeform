#include "zxtimestepper.h"

zxTimeStepper::zxTimeStepper()
{

}

zxTimeStepper::zxTimeStepper(zxForceModel::Ptr forcemodel)
{
    m_forceModel = forcemodel;
    m_dim = m_forceModel->getForceDimension();
    m_residual = Eigen::VectorXd(m_dim);
    m_residual.setZero();
    m_ext_force = m_int_force = m_pos  = m_pos0 = m_vel = m_vel0 = m_dpos = m_residual;
}

void zxTimeStepper::update_mesh()
{
    zxMesh::Ptr mesh = m_forceModel->get_mesh();
    for(size_t i = 0; i < mesh->get_num_nodes(); i++)
    {
        zxNode::Ptr node = mesh->get_node(i);
        node->rt = node->r0 + vec3d(m_pos[3 * i + 0],m_pos[3 * i + 1],m_pos[3 * i + 2]);
    }
}
