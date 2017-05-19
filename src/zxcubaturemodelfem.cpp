#include "zxcubaturemodelfem.h"

zxCubatureModelFEM::zxCubatureModelFEM(zxForceModel::Ptr forcemodel,Eigen::MatrixXd& U,Eigen::VectorXd& Lambda)
{
    m_force_model = forcemodel;

    m_mode_U = U;
    mLamda = Lambda;
    m_mode_UTrans = U.transpose();

    m_dimX = m_dimY = U.cols();
}

Eigen::VectorXd     zxCubatureModelFEM::evaluateY(Eigen::VectorXd& X)
{
    Eigen::VectorXd u = m_mode_U * X;
    Eigen::VectorXd force(u.rows());

    m_force_model->updatePosition(u);
    m_force_model->computeForce(force);

    return m_mode_UTrans * force;
}

Eigen::VectorXd     zxCubatureModelFEM::get_element_column(int el)
{
    Eigen::VectorXd forceVec(get_Y_dim() * get_num_samples());

    Eigen::Map<Eigen::MatrixXd> force(forceVec.data(),get_Y_dim() , get_num_samples());
    force.setZero();
    std::vector<size_t> node_id = m_force_model->get_element_node_id(el);
    for(int isample = 0; isample < get_num_samples(); isample++)
    {
        Eigen::VectorXd ue(node_id.size() * 3), fe(node_id.size() * 3);

        ue.setZero(); fe.setZero();
        for(size_t ii = 0; ii < node_id.size() ; ii++)
            for(int jj = 0; jj < 3; jj++)
                for(int ir = 0; ir < get_X_dim(); ir++)
                ue[3 * ii + jj] += mSampleX(ir,isample) * m_mode_U(3 * node_id[ii] + jj,ir);
        m_force_model->computeElementForce(el,fe,ue);

        for(size_t ii = 0; ii < node_id.size() ; ii++)
            for(int jj = 0; jj < 3; jj++)
                for(int ir = 0; ir < get_Y_dim(); ir++)
                 force(ir,isample) += fe[3 * ii + jj] * m_mode_U(3 * node_id[ii] + jj,ir) * mSampleMagnitudeScale[isample];

    }

    return forceVec;

}

size_t              zxCubatureModelFEM::get_num_elements()
{
    return m_force_model->get_num_elements();
}
