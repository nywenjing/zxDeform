#ifndef ZXAABB_H
#define ZXAABB_H

#include "zxsettings.h"

class zxAABB
{
    ZX_MAKE_SHARED_MACO(zxAABB)
public:
    zxAABB(){reset();}

    static Ptr create() {return Ptr (new zxAABB());}

public:
    bool overlap(Ptr box2,real margin = 0);
    void reset(){ m_lowerB = vec3d(1e60,1e60,1e60); m_upperB = -m_lowerB;}
    void combine(zxAABB::Ptr box2);
    void combine(const vec3d& r);
    size_t get_largest_span_axis();
    real   get_span(size_t axis){return m_upperB[axis] - m_lowerB[axis];}
    const vec3d&  get_lowerB() const {return m_lowerB;}
    const vec3d&  get_upperB() const {return m_upperB;}
    vec3d  get_center(){return 0.5 * (m_lowerB + m_upperB);}
protected:
    vec3d  m_lowerB,m_upperB;
};


#endif // ZKAABB_H
