#ifndef ZXMESH_H
#define ZXMESH_H

#include "zxbasicgeometry.h"
#include "zxabqreader.h"
class zxMesh
{
    ZX_MAKE_SHARED_MACO(zxMesh)
public:
    zxMesh(){}
public:
    size_t                      get_num_elements(){return m_elements.size();}
    size_t                      get_num_nodes(){return m_nodes.size();}
    zxElement::Ptr              get_element(size_t id){return m_elements.at(id);}
    zxElement::Ptr              get_element(size_t id) const{ return m_elements.at(id);}

    zxNode::Ptr                 get_node(size_t id){return m_nodes.at(id);}
    zxNode::ConstPtr            get_node(size_t id) const{return m_nodes.at(id);}

    zxNode::Ptr                 get_node(size_t e_id,size_t id){return get_element(e_id)->get_node(id);}
    zxNode::ConstPtr            get_node(size_t e_id,size_t id) const{return get_element(e_id)->get_node(id);}

    std::vector<zxNode::Ptr>&   get_nodes(){return m_nodes;}
    std::vector<zxElement::Ptr>&get_elements(){return m_elements;}

public:
    std::vector<zxElement::Ptr> m_elements;
    std::vector<zxNode::Ptr>     m_nodes;
};


class zxSolidMesh : public zxMesh
{
ZX_MAKE_SHARED_MACO(zxSolidMesh)
public:
    virtual void buildSurface(Eigen::MatrixXd& V,Eigen::MatrixXi& W_id,Eigen::MatrixXd& W_val,Eigen::MatrixXi& F,int tsLevel = 0) {}

};

class zxSurfaceMesh: public zxMesh
{

};

class zxTetrahedralMesh : public zxSolidMesh
{
ZX_MAKE_SHARED_MACO(zxTetrahedralMesh)
public:
    zxTetrahedralMesh(){}
    zxTetrahedralMesh(std::string filename);

    void    convertToC3D10();
    virtual void buildSurface(Eigen::MatrixXd& V,Eigen::MatrixXi& W_id,Eigen::MatrixXd& W_val,Eigen::MatrixXi& F,int tsLevel);

    static Ptr create(std::string filename) {return std::make_shared<zxTetrahedralMesh>(filename);}
};

class zxHexahedralMesh : public zxSolidMesh
{
    ZX_MAKE_SHARED_MACO(zxHexahedralMesh)
    public:
        zxHexahedralMesh(std::string filename);

        virtual void buildSurface(Eigen::MatrixXd& V,Eigen::MatrixXi& W_id,Eigen::MatrixXd& W_val,Eigen::MatrixXi& F,int tsLevel);

        static Ptr create(std::string filename) {return (Ptr) new zxHexahedralMesh(filename);}

};

class zxTriangularMesh : public zxSurfaceMesh
{
    ZX_MAKE_SHARED_MACO(zxTriangularMesh)

};

#endif // ZXMESH_H
