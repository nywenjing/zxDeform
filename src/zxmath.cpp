#include "zxmath.h"
#include <QMatrix4x4>
void zxModified_SVD(const mat3d& M,mat3d& U,vec3d& Lamda,mat3d& V)
{
    Eigen::JacobiSVD<mat3d> svd(M,Eigen::ComputeFullU | Eigen::ComputeFullV);
    Lamda = svd.singularValues();
    U = svd.matrixU();
    V = svd.matrixV();

    if(U.determinant() < 0)
    {
        if(V.determinant() < 0)
        {
            U *= -1;
            V *= -1;
            return;
        }
        size_t axis = 0;
        for(size_t i = 1; i < 3; i++)
            if(Lamda[axis] < Lamda(i))
                axis = i;
        Lamda[axis] *= -1;
        for(size_t i = 0; i < 3; i++)
            U(i,axis) *= -1;
    }
}

void zx_vega_ComputeDiagonalPstress(zxMaterial::ConstPtr material, const vec3d& lambda,vec3d& diag_P)
{
    vec3d invariants;

    real lambda2[3] = { lambda[0] * lambda[0], lambda[1] * lambda[1], lambda[2] * lambda[2] };
    real IC = lambda2[0] + lambda2[1] + lambda2[2];
    real IIC = lambda2[0] * lambda2[0] + lambda2[1] * lambda2[1] + lambda2[2] * lambda2[2];
    real IIIC = lambda2[0] * lambda2[1] * lambda2[2];

    invariants[0] = IC;
    invariants[1] = IIC;
    invariants[2] = IIIC;

    vec3d dPsidI;

    material->compute_gradient(invariants, dPsidI);

    // PDiag = [ dI / dlambda ]^T * dPsidI

    real mat[9];
    mat[0] = 2.0 * lambda[0];
    mat[1] = 2.0 * lambda[1];
    mat[2] = 2.0 * lambda[2];
    mat[3] = 4.0 * lambda[0] * lambda[0] * lambda[0];
    mat[4] = 4.0 * lambda[1] * lambda[1] * lambda[1];
    mat[5] = 4.0 * lambda[2] * lambda[2] * lambda[2];
    mat[6] = 2.0 * lambda[0] * lambda2[1] * lambda2[2];
    mat[7] = 2.0 * lambda[1] * lambda2[0] * lambda2[2];
    mat[8] = 2.0 * lambda[2] * lambda2[0] * lambda2[1];

    mat3d matM;

    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            matM(i,j) = mat[3 * i + j];
    diag_P = matM.transpose() * dPsidI;
}

int teranToRowMajorMatrix[9] = {0,4,8,1,3,2,6,5,7};
int rowMajorMatrixToTeran[9] = {0,3,5,4,1,7,6,8,2};

int tensor9x9Index(int i, int j, int m, int n)
{
    /*
      |  dP_0/dF_0  dP_0/dF_4  dP_0/dF_8  ...  dP_0/dF_5  |
      |  dP_4/dF_0  dP_4/dF_4  dP_4/dF_8  ...  dP_4/dF_5  |
      |                         ...                       |
      |  dP_5/dF_0  dP_5/dF_4  dP_5/dF_8  ...  dP_5/dF_5  |
     */
    int rowIndex_in9x9Matrix = rowMajorMatrixToTeran[3 * i + j];
    int columnIndex_in9x9Matrix = rowMajorMatrixToTeran[3 * m + n];
    return (9 * rowIndex_in9x9Matrix + columnIndex_in9x9Matrix);
}


double gammaValue(int i, int j, vec3d& sigma, vec3d& invariants, vec3d& gradient, double hessian[6])
{
    /*
     The hessian is in order (11,12,13,22,23,33)
     | 11 12 13 |   | 0 1 2 |
     | 21 22 23 | = | 1 3 4 |
     | 31 32 33 |   | 2 4 5 |
   */

    double tempGammaVec1[3];
    tempGammaVec1[0] = 2.0 * sigma[i];
    tempGammaVec1[1] = 4.0 * sigma[i] * sigma[i] * sigma[i];
    tempGammaVec1[2] = 2.0 * invariants[2] / sigma[i];
    double tempGammaVec2[3];
    tempGammaVec2[0] = 2.0 * sigma[j];
    tempGammaVec2[1] = 4.0 * sigma[j] * sigma[j] * sigma[j];
    tempGammaVec2[2] = 2.0 * invariants[2] / sigma[j];
    double productResult[3];
    productResult[0] = (tempGammaVec2[0] * hessian[0] + tempGammaVec2[1] * hessian[1] +
            tempGammaVec2[2] * hessian[2]);
    productResult[1] = (tempGammaVec2[0] * hessian[1] + tempGammaVec2[1] * hessian[3] +
            tempGammaVec2[2] * hessian[4]);
    productResult[2] = (tempGammaVec2[0] * hessian[2] + tempGammaVec2[1] * hessian[4] +
            tempGammaVec2[2] * hessian[5]);
    return (tempGammaVec1[0] * productResult[0] + tempGammaVec1[1] * productResult[1] +
            tempGammaVec1[2] * productResult[2] + 4.0 * invariants[2] * gradient[2] / (sigma[i] * sigma[j]));
}

void FixPositiveIndefiniteness(double & A11, double & A12, double & A13, double & A22, double & A23, double & A33)
{
    mat3d mat;

    mat(0,0) = A11; mat(0,1) = A12; mat(0,2) = A13;
    mat(1,0) = A12; mat(0,1) = A22; mat(0,2) = A23;
    mat(2,0) = A13; mat(0,1) = A23; mat(0,2) = A33;

    vec3d eigenValues;
    vec3d eigenVectors[3];
    //eigen_sym(mat, eigenValues, eigenVectors);

    Eigen::EigenSolver<mat3d> eigensolver;
    eigensolver.compute(mat);
    eigenValues = eigensolver.eigenvalues().real();
    mat3d eigenvec = eigensolver.eigenvectors().real();

    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
            eigenVectors[i][j] = eigenvec(j,i);

    bool hasNegativeEigenValues = false;
    for(int i = 0; i < 3; i++)
    {
        if (eigenValues[i] < 0)
        {
            hasNegativeEigenValues = true;
            eigenValues[i] = 0;
        }
    }

    if (hasNegativeEigenValues)
    {
        mat3d V;

        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
                V(i,j) = eigenVectors[j][i];

        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
                V(i,j) *= eigenValues[j];



        mat3d newMat = V * V.transpose();
        A11 = newMat(0,0);
        A22 = newMat(1,1);
        A33 = newMat(2,2);
        A12 = newMat(0,1);
        A13 = newMat(0,2);
        A23 = newMat(1,2);
    }
}

void FixPositiveIndefiniteness(double & B11, double & B12)
{
    vec2d eigenValues;
    eigenValues[0] = (B11 - B12);
    eigenValues[1] = (B11 + B12);

    bool hasNegativeEigenValues = false;
    for(int i = 0; i < 2; i++)
    {
        if (eigenValues[i] < 0)
            hasNegativeEigenValues = true;
    }

    if (hasNegativeEigenValues)
    {
        if (eigenValues[0] < 0)
            eigenValues[0] = 0;

        if (eigenValues[1] < 0)
            eigenValues[1] = 0;

        B11 = 0.5 * ( eigenValues[0] + eigenValues[1]);
        B12 = 0.5 * (-eigenValues[0] + eigenValues[1]);
    }
}

void zx_vega_ComputePartialPstress_to_DefGrad(zxMaterial::ConstPtr material, const mat3d U,vec3d sigma,mat3d V,Eigen::MatrixXd& dPsdF_diag,int clamped,bool enforceSPD)
{
    double sigma1square = sigma[0] * sigma[0];
    double sigma2square = sigma[1] * sigma[1];
    double sigma3square = sigma[2] * sigma[2];

    vec3d invariants;
    invariants[0] = sigma1square + sigma2square + sigma3square;
    invariants[1] = (sigma1square * sigma1square +
                     sigma2square * sigma2square +
                     sigma3square * sigma3square);
    invariants[2] = sigma1square * sigma2square * sigma3square;

    //double E[3];
    //E[0] = 0.5 * (Fhats[el][0] * Fhats[el][0] - 1);
    //E[1] = 0.5 * (Fhats[el][1] * Fhats[el][1] - 1);
    //E[2] = 0.5 * (Fhats[el][2] * Fhats[el][2] - 1);

    vec3d gradient;
    material->compute_gradient(invariants, gradient);

    /*
      in order (11,12,13,22,23,33)
      | 11 12 13 |   | 0 1 2 |
      | 21 22 23 | = | 1 3 4 |
      | 31 32 33 |   | 2 4 5 |
    */
    double hessian[6];
    material->compute_hessian(invariants, hessian);

    // modify hessian to compute correct values if in the inversion handling regime
    if (clamped & 1) // first lambda was clamped (in inversion handling)
    {
        hessian[0] = hessian[1] = hessian[2] = 0.0;
    }

    if (clamped & 2) // second lambda was clamped (in inversion handling)
    {
        hessian[1] = hessian[3] = hessian[4] = 0.0;
    }

    if (clamped & 4) // third lambda was clamped (in inversion handling)
    {
        hessian[0] = hessian[1] = hessian[2] = hessian[4] = hessian[5] = 0.0;
    }

    double alpha11 = 2.0 * gradient[0] + 8.0 * sigma1square * gradient[1];
    double alpha22 = 2.0 * gradient[0] + 8.0 * sigma2square * gradient[1];
    double alpha33 = 2.0 * gradient[0] + 8.0 * sigma3square * gradient[1];
    double alpha12 = 2.0 * gradient[0] + 4.0 * (sigma1square+sigma2square) * gradient[1];
    double alpha13 = 2.0 * gradient[0] + 4.0 * (sigma1square+sigma3square) * gradient[1];
    double alpha23 = 2.0 * gradient[0] + 4.0 * (sigma2square+sigma3square) * gradient[1];

    double beta11 = 4.0 * sigma1square * gradient[1] - (2.0 * invariants[2] * gradient[2]) / sigma1square;
    double beta22 = 4.0 * sigma2square * gradient[1] - (2.0 * invariants[2] * gradient[2]) / sigma2square;
    double beta33 = 4.0 * sigma3square * gradient[1] - (2.0 * invariants[2] * gradient[2]) / sigma3square;
    double beta12 = 4.0 * sigma[0] * sigma[1] * gradient[1] - (2.0 * invariants[2] * gradient[2]) / (sigma[0] * sigma[1]);
    double beta13 = 4.0 * sigma[0] * sigma[2] * gradient[1] - (2.0 * invariants[2] * gradient[2]) / (sigma[0] * sigma[2]);
    double beta23 = 4.0 * sigma[1] * sigma[2] * gradient[1] - (2.0 * invariants[2] * gradient[2]) / (sigma[1] * sigma[2]);

    double gamma11 = gammaValue(0, 0, sigma, invariants, gradient, hessian);
    double gamma22 = gammaValue(1, 1, sigma, invariants, gradient, hessian);
    double gamma33 = gammaValue(2, 2, sigma, invariants, gradient, hessian);
    double gamma12 = gammaValue(0, 1, sigma, invariants, gradient, hessian);
    double gamma13 = gammaValue(0, 2, sigma, invariants, gradient, hessian);
    double gamma23 = gammaValue(1, 2, sigma, invariants, gradient, hessian);

    double x1111, x2222, x3333;
    double x2211, x3311, x3322;
    double x2121, x3131, x3232;
    double x2112, x3113, x3223;

    x1111 = alpha11 + beta11 + gamma11;
    x2222 = alpha22 + beta22 + gamma22;
    x3333 = alpha33 + beta33 + gamma33;

    x2211 = gamma12;
    x3311 = gamma13;
    x3322 = gamma23;

    x2121 = alpha12;
    x3131 = alpha13;
    x3232 = alpha23;

    x2112 = beta12;
    x3113 = beta13;
    x3223 = beta23;

    if (enforceSPD)
    {
        FixPositiveIndefiniteness(x1111, x2211, x3311, x2222, x3322, x3333);
        FixPositiveIndefiniteness(x2121, x2112);
        FixPositiveIndefiniteness(x3131, x3113);
        FixPositiveIndefiniteness(x3232, x3223);
    }

    double dPdF_atFhat[81];
    for(int i = 0; i < 81; i++)
        dPdF_atFhat[i] = 0;
    dPdF_atFhat[tensor9x9Index(0,0,0,0)] = x1111;
    dPdF_atFhat[tensor9x9Index(0,0,1,1)] = x2211;
    dPdF_atFhat[tensor9x9Index(0,0,2,2)] = x3311;

    dPdF_atFhat[tensor9x9Index(1,1,0,0)] = x2211;
    dPdF_atFhat[tensor9x9Index(1,1,1,1)] = x2222;
    dPdF_atFhat[tensor9x9Index(1,1,2,2)] = x3322;

    dPdF_atFhat[tensor9x9Index(2,2,0,0)] = x3311;
    dPdF_atFhat[tensor9x9Index(2,2,1,1)] = x3322;
    dPdF_atFhat[tensor9x9Index(2,2,2,2)] = x3333;

    dPdF_atFhat[tensor9x9Index(0,1,0,1)] = x2121;
    dPdF_atFhat[tensor9x9Index(0,1,1,0)] = x2112;

    dPdF_atFhat[tensor9x9Index(1,0,0,1)] = x2112;
    dPdF_atFhat[tensor9x9Index(1,0,1,0)] = x2121;

    dPdF_atFhat[tensor9x9Index(0,2,0,2)] = x3131;
    dPdF_atFhat[tensor9x9Index(0,2,2,0)] = x3113;

    dPdF_atFhat[tensor9x9Index(2,0,0,2)] = x3113;
    dPdF_atFhat[tensor9x9Index(2,0,2,0)] = x3131;

    dPdF_atFhat[tensor9x9Index(1,2,1,2)] = x3232;
    dPdF_atFhat[tensor9x9Index(1,2,2,1)] = x3223;

    dPdF_atFhat[tensor9x9Index(2,1,1,2)] = x3223;
    dPdF_atFhat[tensor9x9Index(2,1,2,1)] = x3232;


    mat3d UT = U.transpose();
    mat3d VT = V.transpose();



    double eiejVector[9];

    for(int i = 0; i < 9; i++)
        eiejVector[i] = 0;
    dPsdF_diag.setZero();
    for (int column=0; column<9; column++)
    {
        eiejVector[column] = 1.0;
        mat3d ei_ej;

        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
                ei_ej(i,j) = eiejVector[3 * i + j];
        mat3d ut_eiej_v = UT*ei_ej*(V);

        double ut_eiej_v_TeranVector[9]; //in Teran order
        ut_eiej_v_TeranVector[rowMajorMatrixToTeran[0]] = ut_eiej_v(0,0);
        ut_eiej_v_TeranVector[rowMajorMatrixToTeran[1]] = ut_eiej_v(0,1);
        ut_eiej_v_TeranVector[rowMajorMatrixToTeran[2]] = ut_eiej_v(0,2);
        ut_eiej_v_TeranVector[rowMajorMatrixToTeran[3]] = ut_eiej_v(1,0);
        ut_eiej_v_TeranVector[rowMajorMatrixToTeran[4]] = ut_eiej_v(1,1);
        ut_eiej_v_TeranVector[rowMajorMatrixToTeran[5]] = ut_eiej_v(1,2);
        ut_eiej_v_TeranVector[rowMajorMatrixToTeran[6]] = ut_eiej_v(2,0);
        ut_eiej_v_TeranVector[rowMajorMatrixToTeran[7]] = ut_eiej_v(2,1);
        ut_eiej_v_TeranVector[rowMajorMatrixToTeran[8]] = ut_eiej_v(2,2);
        double dPdF_resultVector[9]; // not in Teran order
        for (int innerRow=0; innerRow<9; innerRow++)
        {
            double tempResult = 0.0;
            for (int innerColumn=0; innerColumn<9; innerColumn++)
            {
                tempResult += dPdF_atFhat[innerRow*9+innerColumn]*
                        ut_eiej_v_TeranVector[innerColumn];
            }
            dPdF_resultVector[teranToRowMajorMatrix[innerRow]] = tempResult;
        }
        mat3d dPdF_resultMatrix;


        for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
                dPdF_resultMatrix(i,j) = dPdF_resultVector[3 * i + j];

        mat3d u_dpdf_vt = (U)*dPdF_resultMatrix*VT;
        dPsdF_diag.data()[column +  0] = u_dpdf_vt(0,0);
        dPsdF_diag.data()[column +  9] = u_dpdf_vt(0,1);
        dPsdF_diag.data()[column + 18] = u_dpdf_vt(0,2);
        dPsdF_diag.data()[column + 27] = u_dpdf_vt(1,0);
        dPsdF_diag.data()[column + 36] = u_dpdf_vt(1,1);
        dPsdF_diag.data()[column + 45] = u_dpdf_vt(1,2);
        dPsdF_diag.data()[column + 54] = u_dpdf_vt(2,0);
        dPsdF_diag.data()[column + 63] = u_dpdf_vt(2,1);
        dPsdF_diag.data()[column + 72] = u_dpdf_vt(2,2);
        // reset
        eiejVector[column] = 0.0;
    }

}



void zx_plane_dir(vec3d n, vec3d &d0, vec3d d1)
{
    size_t axis = 0;
    for(size_t i = 1; i < 3; i++)
        if(n[i] < n[axis])
            axis = i;


    d0.setZero();
    d0[axis] = 1.0;
    d1 = n.cross(d0);
    d1.normalize();
    d0 = n.cross(d1);
    d0.normalize();
}

double g_degen_normal_epsilon = 1e-6;
double g_cubic_solver_tol = 1e-8;
double g_collision_epsilon = 1e-6;

double clamp(const double &x, const double &a, const double &b)
{
    if(x < a)
        return a;
    if(x > b)
        return b;

    return x;
}

void swap(double& a,double& b)
{
    double tem = a;
    a = b;
    b = tem;
}

double triple(const vec3d &a, const vec3d &b, const vec3d &c)
{
    return a[0]*(b[1]*c[2]-b[2]*c[1])
            +a[1]*(b[2]*c[0]-b[0]*c[2])
            +a[2]*(b[0]*c[1]-b[1]*c[0]);
}

double signed_volume(const vec3d &x0, const vec3d &x1, const vec3d &x2, const vec3d &x3)
{
    // Equivalent to triple(x1-x0, x2-x0, x3-x0), six times the signed volume of the tetrahedron.
    // But, for robustness, we want the result (up to sign) to be independent of the ordering.
    // And want it as accurate as possible...
    // But all that stuff is hard, so let's just use the common assumption that all coordinates are >0,
    // and do something reasonably accurate in fp.

    // This formula does almost four times too much multiplication, but if the coordinates are non-negative
    // it suffers in a minimal way from cancellation error.
    return ( x0[0]*(x1[1]*x3[2]+x3[1]*x2[2]+x2[1]*x1[2])
            +x1[0]*(x2[1]*x3[2]+x3[1]*x0[2]+x0[1]*x2[2])
            +x2[0]*(x3[1]*x1[2]+x1[1]*x0[2]+x0[1]*x3[2])
            +x3[0]*(x1[1]*x2[2]+x2[1]*x0[2]+x0[1]*x1[2]) )

            - ( x0[0]*(x2[1]*x3[2]+x3[1]*x1[2]+x1[1]*x2[2])
            +x1[0]*(x3[1]*x2[2]+x2[1]*x0[2]+x0[1]*x3[2])
            +x2[0]*(x1[1]*x3[2]+x3[1]*x0[2]+x0[1]*x1[2])
            +x3[0]*(x2[1]*x1[2]+x1[1]*x0[2]+x0[1]*x2[2]) );
}

template<class T>
void add_unique(std::vector<T>& a, T e)
{
    for(unsigned int i=0; i<a.size(); ++i)
        if(a[i]==e) return;
    a.push_back(e);
}

void zx_check_point_edge_proximity(bool update, const vec3d &x0, const vec3d &x1, const vec3d &x2,
                                                    double &distance, double &s, vec3d &normal, double normal_multiplier);

void zx_check_point_triangle_proximity(const vec3d &x0, const vec3d &x1, const vec3d &x2, const vec3d &x3, double &distance, double &s1, double &s2, double &s3, vec3d &normal)
{
    vec3d x13=x1-x3;
    double r00=(x13).norm()+1e-30;
    x13/=r00;
    vec3d x23=x2-x3;
    double r01=x23.dot(x13);
    x23-=r01*x13;
    double r11=(x23).norm()+1e-30;
    x23/=r11;
    vec3d x03=x0-x3;
    s2=x23.dot(x03)/r11;
    s1=(x13.dot(x03)-r01*s2)/r00;
    s3=1-s1-s2;
    // check if we are in range
    if(s1>=0 && s2>=0 && s3>=0){
        normal=x0-(s1*x1+s2*x2+s3*x3);
        distance=(normal).norm();
        if(distance>0) normal/=distance;
        else{
            normal=x2-x1.cross(x3-x1);
            normal/=(normal).norm()+1e-30;
        }
    }else{
        double s, d;
        if(s1>0){ // rules out edge 2-3
            zx_check_point_edge_proximity(false, x0, x1, x2, distance, s, normal, 1.);
            s1=s; s2=1-s; s3=0; d=distance;
            zx_check_point_edge_proximity(true, x0, x1, x3, distance, s, normal, 1.);
            if(distance<d){
                s1=s; s2=0; s3=1-s;
            }
        }else if(s2>0){ // rules out edge 1-3
            zx_check_point_edge_proximity(false, x0, x1, x2, distance, s, normal, 1.);
            s1=s; s2=1-s; s3=0; d=distance;
            zx_check_point_edge_proximity(true, x0, x2, x3, distance, s, normal, 1.);
            if(distance<d){
                s1=0; s2=s; s3=1-s; d=distance;
            }
        }else{ // s3>0: rules out edge 1-2
            zx_check_point_edge_proximity(false, x0, x2, x3, distance, s, normal, 1.);
            s1=0; s2=s; s3=1-s; d=distance;
            zx_check_point_edge_proximity(true, x0, x1, x3, distance, s, normal, 1.);
            if(distance<d){
                s1=s; s2=0; s3=1-s;
            }
        }
    }

}


void zx_check_point_edge_proximity(bool update, const vec3d &x0, const vec3d &x1, const vec3d &x2,
                                                    double &distance, double &s, vec3d &normal, double normal_multiplier)
{
    vec3d dx(x2-x1);
    double m2=dx.squaredNorm();
    if(update){
        // find parameter value of closest point on segment
        double this_s=clamp((x2-x0).dot(dx)/m2, 0., 1.);
        // and find the distance
        vec3d this_normal=x0-(this_s*x1+(1-this_s)*x2);
        double this_distance=(this_normal).norm();
        if(this_distance<distance){
            s=this_s;
            distance=this_distance;
            normal=(normal_multiplier/(this_distance+1e-30))*this_normal;
        }
    }else{
        // find parameter value of closest point on segment
        s=clamp((x2-x0).dot(dx)/m2, 0., 1.);
        // and find the distance
        normal=x0-(s*x1+(1-s)*x2);
        distance=(normal).norm();
        normal*=normal_multiplier/(distance+1e-30);
    }
}



void zx_check_edge_edge_proximity(
        const vec3d &x0, const vec3d &x1,
        const vec3d &x2, const vec3d &x3,
        double &distance,
        double &s0, double &s2, vec3d &normal)
{
    // let's do it the QR way for added robustness
    vec3d x01=x0-x1;
    double r00=(x01).norm()+1e-30;
    x01/=r00;
    vec3d x32=x3-x2;
    double r01=x32.dot(x01);
    x32-=r01*x01;
    double r11=(x32).norm()+1e-30;
    x32/=r11;
    vec3d x31=x3-x1;
    s2=x32.dot(x31)/r11;
    s0=(x01.dot(x31)-r01*s2)/r00;
    // check if we're in range
    if(s0<0){
        if(s2<0){
            // check both x1 against 2-3 and 3 against 0-1
            zx_check_point_edge_proximity(false, x1, x2, x3, distance, s2, normal, 1.);
            zx_check_point_edge_proximity(true, x3, x0, x1, distance, s0, normal, -1.);
        }else if(s2>1){
            // check both x1 against 2-3 and 2 against 0-1
            zx_check_point_edge_proximity(false, x1, x2, x3, distance, s2, normal, 1.);
            zx_check_point_edge_proximity(true, x2, x0, x1, distance, s0, normal, -1.);
        }else{
            s0=0;
            // check x1 against 2-3
            zx_check_point_edge_proximity(false, x1, x2, x3, distance, s2, normal, 1.);
        }
    }else if(s0>1){
        if(s2<0){
            // check both x0 against 2-3 and 3 against 0-1
            zx_check_point_edge_proximity(false, x0, x2, x3, distance, s2, normal, 1.);
            zx_check_point_edge_proximity(true, x3, x0, x1, distance, s0, normal, -1.);
        }else if(s2>1){
            // check both x0 against 2-3 and 2 against 0-1
            zx_check_point_edge_proximity(false, x0, x2, x3, distance, s2, normal, 1.);
            zx_check_point_edge_proximity(true, x2, x0, x1, distance, s0, normal, -1.);
        }else{
            s0=1;
            // check x0 against 2-3
            zx_check_point_edge_proximity(false, x0, x2, x3, distance, s2, normal, 1.);
        }
    }else{
        if(s2<0){
            s2=0;
            // check x3 against 0-1
            zx_check_point_edge_proximity(false, x3, x0, x1, distance, s0, normal, -1.);
        }else if(s2>1){
            s2=1;
            // check x2 against 0-1
            zx_check_point_edge_proximity(false, x2, x0, x1, distance, s0, normal, -1.);
        }else{ // we already got the closest points!
            normal=(s0*x0+(1-s0)*x1)-(s2*x2+(1-s2)*x3);
            distance=(normal).norm();
            if(distance>0) normal/=distance;
            else{
                normal=(x1-x0).cross(x3-x2);
                normal/=(normal).norm()+1e-30;
            }
        }
    }

}

void find_coplanarity_times(const vec3d &x0, const vec3d &x1, const vec3d &x2, const vec3d &x3,
                            const vec3d &xnew0, const vec3d &xnew1, const vec3d &xnew2, const vec3d &xnew3,
                            std::vector<double> &possible_times,bool simplex_verbose = false)
{

    if ( simplex_verbose )
    {
        std::cout << "finding coplanarity times... " << std::endl;
    }

    possible_times.clear();

    // cubic coefficients, A*t^3+B*t^2+C*t+D (for t in [0,1])
    vec3d x03=x0-x3, x13=x1-x3, x23=x2-x3;
    vec3d v03=(xnew0-xnew3)-x03, v13=(xnew1-xnew3)-x13, v23=(xnew2-xnew3)-x23;
    double A=triple(v03,v13,v23),
            B=triple(x03,v13,v23)+triple(v03,x13,v23)+triple(v03,v13,x23),
            C=triple(x03,x13,v23)+triple(x03,v13,x23)+triple(v03,x13,x23),
            D=triple(x03,x13,x23);
    const double convergence_tol=g_cubic_solver_tol*(std::fabs(A)+std::fabs(B)+std::fabs(C)+std::fabs(D));

    // find intervals to check, or just solve it if it reduces to a quadratic =============================
    std::vector<double> interval_times;
    double discriminant=B*B-3*A*C; // of derivative of cubic, 3*A*t^2+2*B*t+C, divided by 4 for convenience
    if(discriminant<=0){ // monotone cubic: only one root in [0,1] possible

        if ( simplex_verbose ) { std::cout << "monotone cubic" << std::endl; }

        // so we just
        interval_times.push_back(0);
        interval_times.push_back(1);
    }else{ // positive discriminant, B!=0
        if(A==0){ // the cubic is just a quadratic, B*t^2+C*t+D ========================================
            discriminant=C*C-4*B*D; // of the quadratic
            if(discriminant<=0){
                double t=-C/(2*B);
                if(t>=-g_cubic_solver_tol && t<=1+g_cubic_solver_tol){
                    t=clamp(t, 0., 1.);
                    if(std::fabs(signed_volume((1-t)*x0+t*xnew0, (1-t)*x1+t*xnew1, (1-t)*x2+t*xnew2, (1-t)*x3+t*xnew3))<convergence_tol)
                        possible_times.push_back(t);
                }
            }else{ // two separate real roots
                double t0, t1;
                if(C>0) t0=(-C-std::sqrt(discriminant))/(2*B);
                else    t0=(-C+std::sqrt(discriminant))/(2*B);
                t1=D/(B*t0);
                if(t1<t0) swap(t0,t1);
                if(t0>=-g_cubic_solver_tol && t0<=1+g_cubic_solver_tol) possible_times.push_back(clamp(t0, 0., 1.));
                if(t1>=-g_cubic_solver_tol && t1<=1+g_cubic_solver_tol) add_unique(possible_times, clamp(t1, 0., 1.));
            }

            if ( simplex_verbose )
            {
                std::cout << "A == 0" << std::endl;
                for ( size_t i = 0; i < possible_times.size(); ++i )
                {
                    std::cout << "possible_time: " << possible_times[i] << std::endl;
                }
                std::cout << std::endl;
            }

            return;
        }else{ // cubic is not monotone: divide up [0,1] accordingly =====================================
            double t0, t1;
            if(B>0) t0=(-B-std::sqrt(discriminant))/(3*A);
            else    t0=(-B+std::sqrt(discriminant))/(3*A);
            t1=C/(3*A*t0);
            if(t1<t0) swap(t0,t1);

            if ( simplex_verbose ) { std::cout << "interval times: " << t0 << ", " << t1 << std::endl; }

            interval_times.push_back(0);
            if(t0>0 && t0<1)
                interval_times.push_back(t0);
            if(t1>0 && t1<1)
                interval_times.push_back(t1);
            interval_times.push_back(1);
        }
    }

    if ( simplex_verbose )
    {
        unsigned int n_samples = 20;
        double dt = 1.0 / (double)n_samples;
        double min_val = 1e30;
        for ( unsigned int i = 0; i < n_samples; ++i )
        {
            double sample_t = dt * i;
            double sample_val = signed_volume((1-sample_t)*x0+sample_t*xnew0,
                                              (1-sample_t)*x1+sample_t*xnew1,
                                              (1-sample_t)*x2+sample_t*xnew2,
                                              (1-sample_t)*x3+sample_t*xnew3);

            std::cout << "sample_val: " << sample_val << std::endl;

            min_val = std::min( min_val, fabs(sample_val) );
        }
        std::cout << "min_val: " << min_val << std::endl;
    }


    // look for roots in indicated intervals ==============================================================
    // evaluate coplanarity more accurately at each endpoint of the intervals
    std::vector<double> interval_values(interval_times.size());
    for(size_t i=0; i<interval_times.size(); ++i){
        double t=interval_times[i];
        interval_values[i]=signed_volume((1-t)*x0+t*xnew0, (1-t)*x1+t*xnew1, (1-t)*x2+t*xnew2, (1-t)*x3+t*xnew3);
        if ( simplex_verbose )
        {
            std::cout << "interval time: " << t << ", value: " << interval_values[i] << std::endl;
        }
    }

    if ( simplex_verbose )
    {
        std::cout << "convergence_tol: " << convergence_tol << std::endl;
    }

    // first look for interval endpoints that are close enough to zero, without a sign change
    for(size_t i=0; i<interval_times.size(); ++i){
        if(interval_values[i]==0){
            possible_times.push_back(interval_times[i]);
        }else if(std::fabs(interval_values[i])<convergence_tol){
            if((i==0 || (interval_values[i-1]>=0 && interval_values[i]>=0) || (interval_values[i-1]<=0 && interval_values[i]<=0))
                    &&(i==interval_times.size()-1 || (interval_values[i+1]>=0 && interval_values[i]>=0) || (interval_values[i+1]<=0 && interval_values[i]<=0))){
                possible_times.push_back(interval_times[i]);
            }
        }
    }
    // and then search in intervals with a sign change
    for(size_t i=1; i<interval_times.size(); ++i){
        double tlo=interval_times[i-1], thi=interval_times[i], tmid;
        double vlo=interval_values[i-1], vhi=interval_values[i], vmid;
        if((vlo<0 && vhi>0) || (vlo>0 && vhi<0)){
            // start off with secant approximation (in case the cubic is actually linear)
            double alpha=vhi/(vhi-vlo);
            tmid=alpha*tlo+(1-alpha)*thi;
            int iteration=0;

            if ( simplex_verbose ) { std::cout << "cubic solver tol: " << 1e-2*convergence_tol << std::endl; }

            for(; iteration<50; ++iteration){
                vmid=signed_volume((1-tmid)*x0+tmid*xnew0, (1-tmid)*x1+tmid*xnew1,
                                   (1-tmid)*x2+tmid*xnew2, (1-tmid)*x3+tmid*xnew3);
                if(std::fabs(vmid)<1e-2*convergence_tol) break;
                if((vlo<0 && vmid>0) || (vlo>0 && vmid<0)){ // if sign change between lo and mid
                    thi=tmid;
                    vhi=vmid;
                }else{ // otherwise sign change between hi and mid
                    tlo=tmid;
                    vlo=vmid;
                }
                if(iteration%2) alpha=0.5; // sometimes go with bisection to guarantee we make progress
                else alpha=vhi/(vhi-vlo); // other times go with secant to hopefully get there fast
                tmid=alpha*tlo+(1-alpha)*thi;
            }
            if ( iteration >= 50 && simplex_verbose )
            {
                std::cout << "cubic solve failed" << std::endl;
            }
            possible_times.push_back(tmid);
        }
    }
    sort(possible_times.begin(), possible_times.end());

    if ( simplex_verbose )
    {
        std::cout << "=================" << std::endl;

        for ( size_t i = 0; i < possible_times.size(); ++i )
        {
            std::cout << "possible_time: " << possible_times[i] << std::endl;
        }
        std::cout << std::endl;
    }
}

void degenerate_get_point_triangle_collision_normal(const vec3d &x0, const vec3d &x1, const vec3d &x2, const vec3d &x3,
                                                    double &s1, double &s2, double &s3,
                                                    vec3d& normal )
{

    // try triangle normal at start
    normal=(x2-x1).cross(x3-x1);
    double m=(normal).norm();
    if(m>(g_degen_normal_epsilon * g_degen_normal_epsilon))
    {
        normal/=m;
    }
    else
    {
        // if that didn't work, try vector between points at the start

        normal=(s1*x1+s2*x2+s3*x3)-x0;
        m=(normal).norm();
        if(m>g_degen_normal_epsilon)
        {
            normal/=m;
        }
        else
        {
            // if that didn't work, boy are we in trouble; just get any non-parallel vector
            vec3d dx=x2-x1;
            if(dx[0]!=0 || dx[1]!=0)
            {
                normal=vec3d(dx[1], -dx[0], 0);

                normal.normalize();
            }
            else
            {
                dx=x3-x1;
                if(dx[0]!=0 || dx[1]!=0)
                {
                    normal=vec3d(dx[1], -dx[0], 0);
                    normal.normalize();
                }
                else
                {
                    normal=vec3d(0.0, 1.0, 0.0); // the last resort
                }
            }
        }
    }
}

bool zx_check_point_triangle_collision(const vec3d &x0, const vec3d &x1, const vec3d &x2, const vec3d &x3,
                                                        const vec3d &xnew0, const vec3d &xnew1, const vec3d &xnew2, const vec3d &xnew3,
                                                        double &s1, double &s2, double &s3, vec3d &normal, double &t, double collision_epsilon)
{
    std::vector<double> possible_times;
    s1 = s2 = s3 = 0;
    t = -1;
    find_coplanarity_times(x0, x1, x2, x3, xnew0, xnew1, xnew2, xnew3, possible_times);

    for(size_t a=0; a<possible_times.size(); ++a){
        t=possible_times[a];
        vec3d xt0=(1-t)*x0+t*xnew0, xt1=(1-t)*x1+t*xnew1, xt2=(1-t)*x2+t*xnew2, xt3=(1-t)*x3+t*xnew3;
        double distance;
        zx_check_point_triangle_proximity(xt0, xt1, xt2, xt3, distance, s1, s2, s3, normal);
        if(distance<collision_epsilon){
            // now figure out a decent normal
            if(distance<1e-2*g_degen_normal_epsilon)
            { // if we don't trust the normal...
                // first try the triangle normal at collision time
                normal=(xt2-xt1).cross(xt3-xt1);
                double m=(normal).norm();
                if(m>(g_degen_normal_epsilon * g_degen_normal_epsilon)){
                    normal/=m;
                }
                else
                {
                    degenerate_get_point_triangle_collision_normal( x0, x1, x2, x3, s1, s2, s3, normal );
                }
            }
            return true;
        }
    }
    return false;
}

void degenerate_get_edge_edge_collision_normal(const vec3d &x0, const vec3d &x1, const vec3d &x2, const vec3d &x3,
                                               double s0, double s2, vec3d& normal )
{

    // if that didn't work, try cross-product of edges at the start
    normal=(x1-x0).cross(x3-x2);
    double m=(normal).norm();
    if(m>(g_degen_normal_epsilon * g_degen_normal_epsilon)){
        normal/=m;
    }else{
        // if that didn't work, try vector between points at the start
        normal=(s2*x2+(1-s2)*x3)-(s0*x0+(1-s0)*x1);
        m=(normal).norm();
        if(m>g_degen_normal_epsilon){
            normal/=m;
        }else{
            // if that didn't work, boy are we in trouble; just get any non-parallel vector
            vec3d dx=x1-x0;
            if(dx[0]!=0 || dx[1]!=0){
                normal=vec3d(dx[1], -dx[0], 0);

                normal.normalize();
            }else{
                dx=x3-x2;
                if(dx[0]!=0 || dx[1]!=0){
                    normal=vec3d(dx[1], -dx[0], 0);
                    normal.normalize();
                }else{
                    normal=vec3d(0.0, 1.0, 0.0); // the last resort
                }
            }
        }
    }

}


bool zx_check_edge_edge_collision(const vec3d &x0, const vec3d &x1, const vec3d &x2, const vec3d &x3,
                                                   const vec3d &xnew0, const vec3d &xnew1, const vec3d &xnew2, const vec3d &xnew3,
                                                   double &s0, double &s2, vec3d &normal, double &t, double collision_epsilon)
{
    std::vector<double> possible_times;
    find_coplanarity_times(x0, x1, x2, x3, xnew0, xnew1, xnew2, xnew3, possible_times);
    for(size_t a=0; a<possible_times.size(); ++a){
        t=possible_times[a];
        vec3d xt0=(1-t)*x0+t*xnew0, xt1=(1-t)*x1+t*xnew1, xt2=(1-t)*x2+t*xnew2, xt3=(1-t)*x3+t*xnew3;
        double distance;
        zx_check_edge_edge_proximity(xt0, xt1, xt2, xt3, distance, s0, s2, normal);
        if(distance<collision_epsilon){
            // now figure out a decent normal
            if(distance<1e-2*g_degen_normal_epsilon){ // if we don't trust the normal...
                // first try the cross-product of edges at collision time
                normal=(xt1-xt0).cross(xt3-xt2);
                double m=(normal).norm();
                if(m>(g_degen_normal_epsilon * g_degen_normal_epsilon)){
                    normal/=m;
                }else
                {
                    degenerate_get_edge_edge_collision_normal( x0, x1, x2, x3, s0, s2, normal );
                }
            }
            return true;
        }
    }
    return false;
}
