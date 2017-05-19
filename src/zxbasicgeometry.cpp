#include "zxbasicgeometry.h"



real zxTriangle::get_area()
{
    zxNode::Ptr n0 = get_node(0);
    zxNode::Ptr n1 = get_node(1);
    zxNode::Ptr n2 = get_node(2);

    vec3d avec = (n0->rt - n1->rt).cross(n0->rt - n2->rt);
    return 0.5 * avec.norm();
}

real zxTetrahedron::get_volume()
{
    zxNode::Ptr n0 = get_node(0);
    zxNode::Ptr n1 = get_node(1);
    zxNode::Ptr n2 = get_node(2);
    zxNode::Ptr n3 = get_node(3);


    return 1.0 / 6 * std::abs( (n0->rt - n3->rt).dot( (n1->rt - n3->rt).cross(n2->rt - n3->rt)));

}

real zxSolidElement::get_volume()
{
    real vol = 0.0;
    for(size_t g_id = 0; g_id < get_num_gaussian_points(); g_id++)
    {
        vol += get_material_point(g_id)->m_detJac0 * get_gaussian_weight(g_id);
    }

    return vol;
}
