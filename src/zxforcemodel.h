#ifndef ZXFORCEMODEL_H
#define ZXFORCEMODEL_H

#include "zxsettings.h"
#include "zxsparsematrix.h"
#include "zxmesh.h"
#include "zxmaterial.h"

class zxForceModel
{
public:
    typedef std::shared_ptr<zxForceModel> Ptr;
public:
    zxForceModel();

public:
    virtual size_t      getForceDimension() {return m_dim;}
    virtual size_t      get_num_elements() = 0;
    virtual bool        isSparse(){return m_isSparse;}
    const zxSparseMatrix&  getStiffMatrixTopology(){ assert(m_isSparse);return m_spK_topology;}

public:
    virtual void        updatePosition(const Eigen::VectorXd& u) = 0;
    virtual void        computeForce(Eigen::VectorXd& force) = 0;
    virtual void        computeTangent(zxSparseMatrix& tangentK) = 0;
    virtual void        computeElementForce(int el, Eigen::VectorXd& force,Eigen::VectorXd& u) = 0;


public:
    virtual zxMesh::Ptr         get_mesh(){return m_mesh;}
    virtual zxMaterial::Ptr     get_material(){return m_material;}
    virtual std::vector<size_t>         get_element_node_id(int el) = 0;
    virtual void        computeTangent(Eigen::MatrixXd& ) = 0;

public:
    virtual bool        is_reduced(){return m_is_reduced;}
    virtual void        set_reduced_model(Eigen::MatrixXd& U,std::vector<std::pair<int,real>>& cubatures);
    virtual Eigen::MatrixXd&    get_subspace_modes(){return m_mode;}
    virtual Eigen::MatrixXd&    get_subspace_trans_modes(){return m_modeT;}


protected:
    Eigen::VectorXd     m_rest_positions;
    Eigen::VectorXd     m_current_positions;
    size_t              m_dim;
    bool                m_isSparse;

protected:
    zxSparseMatrix      m_spK_topology;

protected:
    bool                                m_is_reduced;
    Eigen::MatrixXd                     m_dsK_topology;
    std::vector<std::pair<int,real>>    m_cubatures;
    Eigen::MatrixXd                     m_mode,m_modeT;
    std::vector<Eigen::MatrixXd>        m_Ue,m_UeT;

protected:
    zxMesh::Ptr         m_mesh;
    zxMaterial::Ptr     m_material;
};


#endif // ZXFORCEMODEL_H
