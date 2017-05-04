#include "zxcontactconstraint.h"

void zxContactConstraint::setcontact(zxContactPoint::Ptr cp)
{
    cp->m_normal *= 1.0;
}

zxContactConstraint::zxContactConstraint(zxContactPoint* cp)
{
    for(size_t i = 0; i < cp->m_verts.size(); i++)
    {
        zxCollisionMesh::Vert::Ptr vert = cp->m_verts[i];
        if(vert->m_nodes.size() != 0)
        {
            for(size_t j = 0; j < vert->m_nodes.size(); j++)
            {
                zxNode::Ptr node = vert->m_nodes[j];
                m_nodes.push_back(node);
                m_jac.push_back(cp->m_normal * cp->m_weights[i] * vert->m_weights[j]);
                m_weights.push_back(cp->m_weights[i] * vert->m_weights[j]);
            }
        }
        else
        {
            zxNode::Ptr node = zxNode::create();
            node->reset_iac();
            node->r0 = vert->x0;
            node->rp = vert->xp;
            node->rt = vert->x;
            m_nodes.push_back(node);
            m_jac.push_back(cp->m_normal * cp->m_weights[i]);
            m_weights.push_back(cp->m_weights[i]);
        }
    }

    m_cb = cp->m_thickness;
    m_upper_limits = zxInfinity;
    m_lower_limits = 0;
    m_lambda = 0;
}
