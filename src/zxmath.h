#ifndef ZXMATH_H
#define ZXMATH_H

#include "zxsettings.h"
#include "zxmaterial.h"

void zxModified_SVD(const mat3d& M,mat3d& U,vec3d& Lamda,mat3d& V);
void zx_vega_ComputeDiagonalPstress(zxMaterial::ConstPtr material, const vec3d& diag_F,vec3d& diag_P);
//void zx_vega_ComputePartialPstress_to_DefGrad(zxMaterial::ConstPtr material, const mat3d U,const vec3d diag_F,const mat3d V,Eigen::MatrixXd dPsdF_diag);
void zx_vega_ComputePartialPstress_to_DefGrad(zxMaterial::ConstPtr material, const mat3d U,vec3d sigma,mat3d V,Eigen::MatrixXd& dPsdF_diag,int clamped = false,bool enforceSPD = false);

mat4d zx_gl_lookAt(const vec3d& eye,const vec3d& cen,const vec3d& up);
mat4d zx_gl_perspective(real fov,real aspect,real znear,real zfar);

void zx_plane_dir(vec3d n,vec3d& d0,vec3d d1);


void zx_check_point_triangle_proximity(const vec3d &x0, const vec3d &x1, const vec3d &x2, const vec3d &x3,
                                    double &distance, double &s1, double &s2, double &s3, vec3d &normal);

void zx_check_edge_edge_proximity(const vec3d &x0, const vec3d &x1, const vec3d &x2, const vec3d &x3,
                               double &distance, double &s0, double &s2, vec3d &normal);

bool zx_check_point_triangle_collision(const vec3d &x0, const vec3d &x1, const vec3d &x2, const vec3d &x3,
                                    const vec3d &xnew0, const vec3d &xnew1, const vec3d &xnew2, const vec3d &xnew3,
                                    double &s1, double &s2, double &s3, vec3d &normal, double &t, double collision_epsilon);

bool zx_check_edge_edge_collision(const vec3d &x0, const vec3d &x1, const vec3d &x2, const vec3d &x3,
                               const vec3d &xnew0, const vec3d &xnew1, const vec3d &xnew2, const vec3d &xnew3,
                               double &s0, double &s2, vec3d &normal, double &t, double collision_epsilon);

#endif // ZXMATH_H
