#include "zxrendermesh.h"

zxRenderMesh::zxRenderMesh()
{
    m_ambient.setZero();
    m_specular.setZero();
    m_diffuse = vec4f(0.8,0.8,0.8,1.0);

}

zxRenderMesh::zxRenderMesh(std::string filename)
    :zxCollisionMesh(filename)
{

    m_ambient.setZero();
    m_specular.setZero();
    m_diffuse = vec4f(0.8,0.8,0.8,1.0);
}

zxRenderMesh::zxRenderMesh(zxSolidMesh::Ptr mesh,size_t tsLevel)
    :zxCollisionMesh(mesh,tsLevel)
{

    m_ambient.setZero();
    m_specular.setZero();
    m_diffuse = vec4f(0.8,0.8,0.8,1.0);
}
