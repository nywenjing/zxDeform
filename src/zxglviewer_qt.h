#ifndef ZXGLVIEWER_QT_H
#define ZXGLVIEWER_QT_H

#include "zxglviewer.h"
#include <QOpenGLWidget>
#include <QMouseEvent>
#include <QTimer>
#include <QMatrix4x4>

class zxGLViewer_QT : public QOpenGLWidget, public zxGLViewer
{
    Q_OBJECT
public:
    zxGLViewer_QT();

public:
    virtual void initGL(){}
    virtual void render(){}
    virtual void init() = 0;
public:
    virtual void initializeGL();
    virtual void paintGL();
    virtual void resizeGL(int w, int h){m_width = w; m_height = h;}

protected:
    virtual void mousePressEvent(QMouseEvent* e);
    virtual void mouseMoveEvent(QMouseEvent* e);
    virtual void wheelEvent(QWheelEvent* e);

public:
    virtual void pressKeyEvent(QKeyEvent* e);
    virtual void toggleAnimated();

public:
    virtual mat4d get_modelview_matrix();
    virtual mat4d get_perspective_matrix();

protected slots:
    virtual void animateSlot(){animate(); update();}

protected:
    QPoint  m_old_mouse_point;
    Qt::MouseButton m_pressed_button;
    QTimer      m_timer;

};

#endif // ZXGLVIEWER_QT_H
