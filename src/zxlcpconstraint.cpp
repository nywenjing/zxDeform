#include "zxlcpconstraint.h"

zxLCPConstraint::zxLCPConstraint()
{
    m_lambda = 0;

}

void zxLCPConstraint::init()
{
    m_invJtAJ = 0;
    for(size_t i = 0; i < m_nodes.size(); i++)
    {
        zxNode::Ptr node = m_nodes[i];
        m_invJtAJ += m_jac[i].dot(node->m_lcp_invA * m_jac[i]);
    }

    m_invJtAJ = 1.0/m_invJtAJ;
}

real zxLCPConstraint::get_residual()
{
    real gap = 0;
    for(size_t i = 0; i < m_nodes.size(); i++)
    {
        zxNode::Ptr node = m_nodes[i];
        gap += m_jac[i].dot(node->rt + node->m_lcp_dx + node->m_lcp_pcg_dx);
    }
    return m_cb - gap;
}

void zxLCPConstraint::reset_projected_gs()
{
    for(size_t i = 0; i < m_nodes.size(); i++)
    {
        zxNode::Ptr node = m_nodes[i];
        node->m_lcp_pcg_dx.setZero();
    }

}

void zxLCPConstraint::do_projected_gs_solve()
{
    real cb = get_residual();

    real rb = 0;

    for(size_t j = 0; j < m_nodes.size(); j++)
    {
        zxNode::Ptr n = m_nodes[j];
        rb += m_jac[j].dot(n->m_lcp_invA * n->m_lcp_res_imp);
    }

    cb -= rb;

    real clamda = m_lambda + m_invJtAJ * cb;
    clamda = std::max(clamda,m_lower_limits);
    clamda = std::min(clamda,m_upper_limits);

    for(size_t i = 0; i < m_nodes.size(); i++)
    {
        zxNode::Ptr node = m_nodes[i];
        node->m_lcp_pcg_dx += node->m_lcp_invA * m_jac[i] * (clamda - m_lambda);
    }

    m_lambda = clamda;

}

void zxLCPConstraint::zero_lcp_impulse()
{
    for(size_t i = 0; i < m_nodes.size(); i++)
    {
        zxNode::Ptr node = m_nodes[i];
        node->m_lcp_imp.setZero();
    }


}

void zxLCPConstraint::set_lcp_impulse()
{
    for(size_t i = 0; i < m_nodes.size(); i++)
    {
        zxNode::Ptr node = m_nodes[i];
        node->m_lcp_imp += m_lambda * m_jac[i];
    }
}
