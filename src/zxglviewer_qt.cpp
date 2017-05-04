#include "zxglviewer_qt.h"

zxGLViewer_QT::zxGLViewer_QT()
{
    m_timer.setParent(this);
    connect(&m_timer, SIGNAL(timeout()), this, SLOT(animateSlot()));
    m_timer.start(0);
    m_timer.stop();


}

void zxGLViewer_QT::initializeGL()
{
    init();

}

void zxGLViewer_QT::paintGL()
{
    render();
}



void zxGLViewer_QT::mousePressEvent(QMouseEvent* e)
{
    m_old_mouse_point = e->pos();
    m_pressed_button = e->button();

}

void zxGLViewer_QT::mouseMoveEvent(QMouseEvent* e)
{
    QPoint e_pos = e->pos();

    if(m_pressed_button == Qt::LeftButton)
    {
        this->rotate(vec2d(m_old_mouse_point.x(),m_old_mouse_point.y()),vec2d(e_pos.x(),e_pos.y()));

    }

    m_old_mouse_point = e_pos;
    update();
}

void zxGLViewer_QT::wheelEvent(QWheelEvent* e)
{
    real s = 0.1;
    if(e->delta() < 0)
        s *= -1;
    zoom(s);


    update();
}

void zxGLViewer_QT::pressKeyEvent(QKeyEvent *e)
{
    if(e->key() == Qt::Key_Space)
        toggleAnimated();
}
void zxGLViewer_QT::toggleAnimated()
{

    zxGLViewer::toggleAnimated();
    if(m_isAnimated)
        m_timer.start();
    else
        m_timer.stop();

    bool flag = (m_isAnimated == true);

    std::cout<<"toggle animation: "<<flag<<std::endl;
}


mat4d zx_gl_lookAt(const vec3d& eye,const vec3d& cen,const vec3d& up)
{


}



mat4d zxGLViewer_QT::get_modelview_matrix()
{
    vec3d eye = m_cen - m_rad * m_forward;
    vec3d up = m_up;
    vec3d cen = m_cen;
    QMatrix4x4 q_mat;
    QVector3D q_eye,q_cen,q_up;
    for(int i = 0; i < 3; i++)
    {
        q_eye[i] = eye[i];
        q_cen[i] = cen[i];
        q_up[i] = up[i];
    }
    q_mat.lookAt(q_eye,q_cen,q_up);

    mat4d mat;
    for(int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++)
            mat(i,j) = q_mat(i,j);

    return mat;

}

mat4d zxGLViewer_QT::get_perspective_matrix()
{
    real fov = 60;
    real aspect = width() * 1.0/ height();
    real znear = 0.1;
    real zfar = 1000;
    QMatrix4x4 q_mat;
    q_mat.perspective(fov,aspect,znear,zfar);

    mat4d mat;
    for(int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++)
            mat(i,j) = q_mat(i,j);

    return mat;

}
