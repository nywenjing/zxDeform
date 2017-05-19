#ifndef ZXSIMVIEWER_H
#define ZXSIMVIEWER_H

#include "zxglviewer_qt.h"
#include "zxsimworld.h"
#include "zxalmsimulator.h"

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
    virtual void init_femsim1();
    virtual void init_femsim2();
    virtual void init_femsim3();



protected:
    zxSimWorld::Ptr m_world;
};

#endif // ZXSIMVIEWER_H
