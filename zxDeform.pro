#-------------------------------------------------
#
# Project created by QtCreator 2017-04-28T14:55:53
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = zxDeform
TEMPLATE = app

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp

INCLUDEPATH += ./src\
./third_party/eigen3.3.3 \
./third_party


SOURCES += main.cpp\
        mainwindow.cpp \
    src/zxmath.cpp \
    src/zxsettings.cpp \
    src/zxbasicgeometry.cpp \
    src/zx_test_program.cpp \
    src/zxmesh.cpp \
    src/zxforcemodel.cpp \
    src/zxsparsematrix.cpp \
    src/zxmaterial.cpp \
    src/zxgaussiantraits.cpp \
    src/zxnonlinearfem_forcemodel_sparse.cpp \
    src/zxtimestepper.cpp \
    src/zxtimestepperbackwardeuler.cpp \
    src/zxaabb.cpp \
    src/zxbody.cpp \
    src/zxsimworld.cpp \
    src/zxcollisionmesh.cpp \
    src/zxbvh.cpp \
    src/zxrendermesh.cpp \
    src/zxglviewer.cpp \
    src/zxglviewer_qt.cpp \
    zxsimviewer.cpp \
    src/zxneohookeanmaterial.cpp \
    src/zxcollisionalgorithm.cpp \
    src/zxcontactpoint.cpp \
    src/zxlcpconstraint.cpp \
    src/zxcontactconstraint.cpp \
    src/zxfrictionconstraint.cpp

HEADERS  += mainwindow.h \
    src/zxmath.h \
    src/zxsettings.h \
    src/zxbasicgeometry.h \
    src/zxmesh.h \
    src/zxforcemodel.h \
    src/zxsparsematrix.h \
    src/zxmaterial.h \
    src/zxgaussiantraits.h \
    src/zxnonlinearfem_forcemodel_sparse.h \
    src/zxtimestepper.h \
    src/zxtimestepperbackwardeuler.h \
    src/zxaabb.h \
    src/zxbody.h \
    src/zxsimworld.h \
    src/zxcollisionmesh.h \
    src/zxbvh.h \
    src/zxrendermesh.h \
    src/zxglviewer.h \
    src/zxglviewer_qt.h \
    zxsimviewer.h \
    src/zxneohookeanmaterial.h \
    src/zxcollisionalgorithm.h \
    src/zxcontactpoint.h \
    src/zxlcpconstraint.h \
    src/zxcontactconstraint.h \
    src/zxfrictionconstraint.h

FORMS    += mainwindow.ui
