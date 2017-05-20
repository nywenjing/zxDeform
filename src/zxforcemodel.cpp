#include "zxforcemodel.h"

zxForceModel::zxForceModel()
{
    m_is_reduced = false;
}

void zxForceModel::set_reduced_model(Eigen::MatrixXd &U, std::vector<std::pair<int, real> > &cubatures)
{

    m_is_reduced = true;
    m_mode = U; m_modeT = U.transpose();
    m_cubatures = cubatures;
    m_dim = m_mode.cols();

    m_Ue.resize(cubatures.size());
    m_UeT.resize(cubatures.size());

    for(size_t el = 0; el < cubatures.size(); el++)
    {
        zxElement::Ptr element = m_mesh->get_element(cubatures[el].first);

        Eigen::MatrixXd& Ue = m_Ue[el];
        Ue.resize(element->get_num_nodes() * 3, U.cols());

        for(int i = 0; i < element->get_num_nodes(); i++)
            for(int j = 0; j < 3; j++)
                for(int ir = 0; ir < U.cols(); ir++)
                    Ue(3 * i + j,ir) = U(3 * element->get_node(i)->m_id + j,ir);

        m_UeT[el] = Ue.transpose();

    }

}
