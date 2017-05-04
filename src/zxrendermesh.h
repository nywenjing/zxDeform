#ifndef ZXRENDERMESH_H
#define ZXRENDERMESH_H

#include "zxcollisionmesh.h"

class zxRenderMesh : public zxCollisionMesh
{
    ZX_MAKE_SHARED_MACO(zxRenderMesh)
    using zxCollisionMesh::Vert;
    using zxCollisionMesh::Edge;
    using zxCollisionMesh::Face;
public:
    zxRenderMesh();
    zxRenderMesh(std::string filename);
    zxRenderMesh(zxSolidMesh::Ptr mesh);

    static Ptr create(std::string filename) {return Ptr (new zxRenderMesh(filename));}
    static Ptr create(zxSolidMesh::Ptr mesh){ return Ptr(new zxRenderMesh(mesh));}



public:
    vec4f m_diffuse,m_specular,m_ambient;
    float m_shiness;
};

#endif // ZXRENDERMESH_H
