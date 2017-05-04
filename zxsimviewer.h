#ifndef ZXSIMVIEWER_H
#define ZXSIMVIEWER_H

#include "zxglviewer_qt.h"
#include "zxsimworld.h"

class zxSimViewer : public zxGLViewer_QT
{
    Q_OBJECT
public:
    zxSimViewer();

public:
    virtual void init();
    virtual void render();
    virtual void animate();

public:
    virtual void init_femsim0();



protected:
    zxSimWorld::Ptr m_world;
};

#endif // ZXSIMVIEWER_H
