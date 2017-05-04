#ifndef ZXSPARSEMATRIX_H
#define ZXSPARSEMATRIX_H

#include "zxsettings.h"

class zxSparseMatrix : public Eigen::SparseMatrix<real>
{
public:
    zxSparseMatrix();

    zxSparseMatrix(int rows,int cols) : Eigen::SparseMatrix<double>(rows,cols){}
    zxSparseMatrix(const Eigen::SparseMatrix<double>& mat)
        :Eigen::SparseMatrix<double>(mat){}
    zxSparseMatrix& operator = (const Eigen::SparseMatrix<double>& mat)
    {
        Eigen::SparseMatrix<double>::operator =(mat);

        return *this;
    }
};

#endif // ZXSPARSEMATRIX_H
