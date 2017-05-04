#ifndef ZXGLVIEWER_H
#define ZXGLVIEWER_H

#include "zxsettings.h"

class zxGLViewer
{
public:
    zxGLViewer();

public:
    virtual void initGL() = 0;
    virtual void render() = 0;
    virtual void init() = 0;
    virtual void animate() {}
public:
    virtual void rotate(vec2d pfrom,vec2d pto);
    virtual void translate(vec2d pfrom,vec2d pto);
    virtual void zoom(real s);
    virtual void showEntireScene();
    virtual void toggleAnimated();
    virtual mat4d get_modelview_matrix() = 0;
    virtual mat4d get_perspective_matrix() = 0;


protected:
    real    m_rad;
    vec3d   m_cen;
    vec3d   m_up;
    vec3d   m_forward;
    vec3d   m_side;
    size_t  m_width;
    size_t  m_height;
    bool    m_isAnimated;

};

#endif // ZXGLVIEWER_H
