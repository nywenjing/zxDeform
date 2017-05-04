#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::keyPressEvent(QKeyEvent *e)
{
    zxGLViewer_QT* viewer = dynamic_cast<zxGLViewer_QT*>(centralWidget());
    if(viewer != nullptr)
    {
        viewer->pressKeyEvent(e);
    }
}
