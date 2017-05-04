#include "zxbvh.h"



zxBVHNode::zxBVHNode(zxAABBData::Ptr element)
{
    m_aabb = zxAABB::create();
    m_aabb->reset();
    m_aabb->combine(element->get_aabb());
    m_data = element;

    m_left = m_right = nullptr;
}

zxBVHNode::zxBVHNode(std::vector<zxAABBData::Ptr> &elements)
{
    m_aabb = zxAABB::create();
    m_aabb->reset();

    for(size_t i = 0; i < elements.size(); i++)
    {
        zxAABBData::Ptr el = elements[i];
        el->update_aabb(false);
        m_aabb->combine(el->get_aabb());
    }

    if(elements.size() == 1)
    {
        m_left = m_right = nullptr;
        m_data = elements[0];
    }
    else if(elements.size() == 2)
    {
        m_left = zxBVHNode::create(elements[0]);
        m_right = zxBVHNode::create(elements[1]);
    }
    else
    {
        std::vector<zxAABBData::Ptr> leftData,rightData;
        size_t axis = m_aabb->get_largest_span_axis();
        real middle = m_aabb->get_center()[axis];

        for(size_t i = 0; i < elements.size(); i++)
        {
            zxAABB::Ptr aabb = elements[i]->get_aabb();
            vec3d cen = aabb->get_center();

            if(cen[axis] < middle)
                leftData.push_back(elements[i]);
            else
                rightData.push_back(elements[i]);
        }

        if(leftData.size() == 0 || rightData.size() == 0)
        {
            leftData.clear();
            rightData.clear();

            size_t halfEl = elements.size() / 2;

            leftData.insert(leftData.begin(),elements.begin(),elements.begin() + halfEl);
            rightData.insert(rightData.begin(),elements.begin() + halfEl,elements.end());
        }

        m_left = zxBVHNode::create(leftData);
        m_right = zxBVHNode::create(rightData);

    }

}

void test_collision(zxBVHNode::Ptr node0,zxBVHNode::Ptr node1)
{

}

void zxBVHNode::collid(zxBVHNode::Ptr other_node, zxBVHCollider::Ptr collider)
{
    real margin = collider->get_margin();
    assert(m_aabb != nullptr);
    if(!m_aabb->overlap(other_node->get_aabb(),collider->get_margin()))
        return;


    if(is_leaf() && other_node->is_leaf())
    {
        collider->process_collision(this,other_node.get());
        assert(m_aabb != nullptr);
    }
    else if(is_leaf())
    {
        collid(other_node->get_left(),collider);
        collid(other_node->get_right(),collider);
    }
    else
    {
        m_left->collid(other_node,collider);
        m_right->collid(other_node,collider);
    }


}

void zxBVHNode::self_collid(zxBVHCollider::Ptr collider)
{
    if(is_leaf())
        return;
    m_left->collid(m_right,collider);
    m_left->self_collid(collider);
    m_right->self_collid(collider);

}

void zxBVHNode::refit()
{
    if(is_leaf())
        m_aabb = m_data->get_aabb();
    else
    {
        m_left->refit();
        m_right->refit();

        m_aabb->reset();
        m_aabb->combine(m_left->get_aabb());
        m_aabb->combine(m_right->get_aabb());
    }
}

zxBVHTree::zxBVHTree(std::vector<zxAABBData::Ptr> &elements)
{
    m_root = zxBVHNode::create(elements);
}
