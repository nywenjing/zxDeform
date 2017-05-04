#include "mainwindow.h"
#include <QApplication>
#include "zxsimviewer.h"

void test_basic_geometry();
int main(int argc, char *argv[])
{
    //test_basic_geometry();

    QApplication a(argc, argv);
    MainWindow w;

    zxSimViewer viewer;
    w.resize(800,600);
    w.setCentralWidget(&viewer);
    w.show();
    return a.exec();
}
