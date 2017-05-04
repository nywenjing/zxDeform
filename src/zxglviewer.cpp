#include "zxglviewer.h"

zxGLViewer::zxGLViewer()
{
    m_cen = vec3d(0.0,0.0,0.0);
    m_rad = 1.0;
    m_up = vec3d(0.0,1.0,0.0);
    m_forward = vec3d(0.0,0.0,-1.0);
    m_side = vec3d(1.0,0.0,0.0);

    m_isAnimated = false;
}

void zxGLViewer::rotate(vec2d pfrom,vec2d pto)
{
    double dx = pto[0] - pfrom[0];
    double dy = pto[1] - pfrom[1];

    vec3d axis = -dx * m_side + dy * m_up;
    axis = axis.cross(m_forward);
    double angle = sqrt(dx * dx + dy * dy) / m_width * M_PI;
    axis.normalize();
    mat3d  rot = Eigen::AngleAxisd(angle,axis).toRotationMatrix();

    m_up = rot * m_up;
    m_side = rot * m_side;
    m_forward = rot * m_forward;
}

void zxGLViewer::translate(vec2d pfrom,vec2d pto)
{
    double dx = pto[0] - pfrom[0];
    double dy = pto[1] - pfrom[1];
    m_cen += (dx * m_side + dy * m_up) * m_rad;
}

void zxGLViewer::zoom(real s)
{
    m_rad -= s * m_rad;
}

void zxGLViewer::showEntireScene()
{

}

void zxGLViewer::toggleAnimated()
{
    m_isAnimated = !m_isAnimated;
}
