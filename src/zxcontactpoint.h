#ifndef ZXCONTACTPOINT_H
#define ZXCONTACTPOINT_H

#include "zxsettings.h"
#include "zxcollisionmesh.h"

class zxContactPoint
{
    ZX_MAKE_SHARED_MACO_DEFAULT_CREAT(zxContactPoint)
public:
    zxContactPoint();

public:
    vec3d m_normal;
    real  m_thickness;

    std::vector<zxCollisionMesh::Vert::Ptr> m_verts;
    std::vector<real>                       m_weights;
};

#endif // ZXCONTACTPOINT_H
