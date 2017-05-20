#include "zxlcpconstraint.h"

zxLCPConstraint::zxLCPConstraint()
{
    m_lambda = 0;
    m_penalty = 1e6;

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

void zxLCPConstraint::set_alm_impulse()
{
    for(size_t i = 0; i < m_nodes.size(); i++)
    {
        zxNode::Ptr node = m_nodes[i];
        real Ln = m_lambda + get_residual() * m_penalty;

        Ln = std::max(Ln,m_lower_limits);
        Ln = std::min(Ln,m_upper_limits);

        node->m_lcp_imp += Ln  * m_jac[i];
    }
}

void zxLCPConstraint::compute_alm_Ab()
{
    real jacb = 0;

    real Ln = m_lambda + get_residual() * m_penalty;

    if(Ln > m_lower_limits && Ln < m_upper_limits)
        Ln = 1.0;
    else
        Ln = 0.0;

    for(size_t i = 0; i < m_nodes.size(); i++)
    {
        zxNode::Ptr node = m_nodes[i];
        jacb += node->m_alm_b.dot(m_jac[i]);

        if(std::isnan(jacb))
            i++;
    }

    jacb *= m_penalty * Ln;

    for(size_t i = 0; i < m_nodes.size(); i++)
    {
        zxNode::Ptr node = m_nodes[i];
        node->m_alm_Ab += jacb * m_jac[i];


    }
}

void zxLCPConstraint::do_alm_augment()
{
    m_lambda += get_residual() * m_penalty;

    m_lambda = std::max(m_lambda,m_lower_limits);
    m_lambda = std::min(m_lambda,m_upper_limits);
}
