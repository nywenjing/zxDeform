#include "zxaabb.h"

bool zxAABB::overlap(Ptr box2,real margin)
{
    for(size_t i = 0; i < 3; i++)
    {
        if(m_lowerB[i] > box2->get_upperB()[i] + margin)
            return false;
        if(m_upperB[i] < box2->get_lowerB()[i] - margin)
            return false;
    }

    return true;

}

void zxAABB::combine(zxAABB::Ptr box2)
{
    m_lowerB = m_lowerB.array().min(box2->get_lowerB().array());
    m_upperB = m_upperB.array().max(box2->get_upperB().array());
}

void zxAABB::combine(const vec3d& r)
{
    m_lowerB = m_lowerB.array().min(r.array());
    m_upperB = m_upperB.array().max(r.array());

}

size_t zxAABB::get_largest_span_axis()
{
    size_t axis = 0;
    real max_span = get_span(0);

    for(size_t i = 1; i < 2; i++)
        if(max_span < get_span(i))
        {
            max_span = get_span(i);
            axis = i;
        }

    return axis;
}
