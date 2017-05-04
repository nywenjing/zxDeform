#include "zxfrictionconstraint.h"

zxFrictionConstraint::zxFrictionConstraint()
{
    m_friction = 0;
}


zxFrictionConstraint::zxFrictionConstraint(zxContactConstraint::Ptr cc,vec3d dir)
{
    m_contact = cc;

    m_nodes = cc->get_nodes();
    m_weights = cc->get_weights();

    for(size_t i = 0; i < m_nodes.size(); i++)
    {
        m_jac.push_back(dir * m_weights[i]);
    }

    m_cb = 0;

    for(size_t i = 0; i < m_nodes.size(); i++)
        m_cb += m_jac[i].dot(m_nodes[i]->rp);
    m_upper_limits = 0;
    m_lower_limits = 0;
    m_lambda = 0;
}

void zxFrictionConstraint::do_projected_gs_solve()
{
    m_upper_limits = m_contact->get_lambda() * m_friction;
    m_lower_limits = -m_upper_limits;
    zxLCPConstraint::do_projected_gs_solve();
}
