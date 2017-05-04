#ifndef ZXBVH_H
#define ZXBVH_H

#include "zxsettings.h"
#include "zxaabb.h"

class zxAABBData
{
    ZX_MAKE_SHARED_MACO_DEFAULT_CREAT(zxAABBData)

    zxAABBData()
    {
        m_aabb = zxAABB::create();
    }

public:
      virtual zxAABB::Ptr get_aabb(){return m_aabb;}
      virtual void update_aabb(bool ccd) {assert(false);}

protected:
    zxAABB::Ptr m_aabb;
};

class zxBVHNode;

class zxBVHCollider
{
    ZX_MAKE_SHARED_MACO(zxBVHCollider)
    zxBVHCollider()
    {
        m_margin = 0;
    }

public:
    virtual void process_collision(zxBVHNode* node0,zxBVHNode* node1) = 0;
    virtual void set_margin(real margin){ m_margin = margin;}
    virtual real get_margin(){return m_margin;}

protected:
    real    m_margin;

};

class zxBVHNode
{
    ZX_MAKE_SHARED_MACO(zxBVHNode)
public:
    zxBVHNode(zxAABBData::Ptr element);
    zxBVHNode(std::vector<zxAABBData::Ptr>& elements);

public:
    static Ptr create(std::vector<zxAABBData::Ptr>& elements){return Ptr (new zxBVHNode(elements)); }
    static Ptr create(zxAABBData::Ptr element){return Ptr (new zxBVHNode(element));}
public:
    void  collid(zxBVHNode::Ptr other_node,zxBVHCollider::Ptr collider);
    void  self_collid(zxBVHCollider::Ptr collider);
    bool  is_leaf(){return m_left == nullptr;}
    zxBVHNode::Ptr get_left(){return m_left;}
    zxBVHNode::Ptr get_right(){return m_right;}
    zxAABB::Ptr get_aabb(){return m_aabb;}
    zxAABBData::Ptr get_data(){return m_data;}

    void              refit();
protected:
    zxAABB::Ptr m_aabb;
    zxBVHNode::Ptr m_left,m_right;
    zxAABBData::Ptr m_data;
};

class zxBVHTree
{
    ZX_MAKE_SHARED_MACO(zxBVHTree)
public:
    zxBVHTree(std::vector<zxAABBData::Ptr>& elements);

    static Ptr create(std::vector<zxAABBData::Ptr>& elements){return Ptr(new zxBVHTree(elements));}
    zxBVHNode::Ptr    get_root(){return m_root;}
    void              refit(){ m_root->refit();}

public:
    zxBVHNode::Ptr  m_root;
};

#endif // ZXBVH_H
