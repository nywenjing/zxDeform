#ifndef SETTINGS_H
#define SETTINGS_H

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"
#include "Eigen/Geometry"
#include "Eigen/Eigenvalues"
#include "Eigen/SVD"
#include <vector>
#include <list>
#include <thread>
#include <numeric>
#include <algorithm>
#include <string>
#include <math.h>
#include <fstream>
#include <iostream>
#include <queue>


typedef double real;
typedef Eigen::Vector3d vec3d;
typedef Eigen::Vector2d vec2d;
typedef Eigen::Vector4d vec4d;
typedef Eigen::Matrix3d mat3d;
typedef Eigen::Matrix4d mat4d;

typedef Eigen::Vector3f vec3f;
typedef Eigen::Vector2f vec2f;
typedef Eigen::Vector4f vec4f;
typedef Eigen::Matrix3f mat3f;
typedef Eigen::Matrix4f mat4f;

#define zxInfinity 1e30
#define zxEPSILON  1e-6

#define ZX_MAKE_SHARED_MACO(T) \
public:\
typedef std::shared_ptr<T> Ptr;\
typedef std::shared_ptr<T const> ConstPtr;

#define ZX_MAKE_SHARED_MACO_DEFAULT_CREAT(T) \
public:\
typedef std::shared_ptr<T> Ptr;\
typedef std::shared_ptr<T const> ConstPtr;\
static Ptr create(){return Ptr(new T());}

#endif // SETTINGS_H
