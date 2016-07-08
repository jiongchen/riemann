#ifndef GEOMETRY_EXTEND_H
#define GEOMETRY_EXTEND_H

#include <zjucad/matrix/matrix.h>
#include <Eigen/Dense>

namespace riemann {

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

void gen_rand_orth_vec_3d(const double *x, double *xT);

int calc_vert_local_frame(const zjucad::matrix::matrix<size_t> &tris,
                          const zjucad::matrix::matrix<double> &nods,
                          zjucad::matrix::matrix<double> &frame);

void project_vector_on_subspace(const double *vect, const size_t dim,
                                const double *basis, const size_t sub_dim,
                                double *proj_vect);

void project_point_on_subspace(const double *x, const size_t dim,
                               const double *origin, const double *basis, const size_t sub_dim,
                               double *proj_x);

void project_point_on_plane(const double *x, const double *origin,
                            const double *tan0, const double *tan1,
                            double *proj_x);

void rotate_and_scale_triangle(const double *x, const double *quaternion, const double s, const double *xnew);

int calc_one_ring_face(const zjucad::matrix::matrix<size_t> &tris,
                       std::vector<std::vector<size_t>> &p2f);

void calc_face_local_frame(const mati_t &tris, const matd_t &nods, matd_t &origin, matd_t &axis);

void calc_local_uv(const mati_t &tris, const matd_t &nods, const matd_t &origin, const matd_t &axis, matd_t &uv);

void calc_tris_cot_value(const mati_t &tris, const matd_t &nods, matd_t &cotv);

inline Eigen::Matrix3d RX(const double alpha) {
  Eigen::Matrix3d R = Eigen::Matrix3d::Zero();
  R(0, 0) = 1;
  R(1, 1) = R(2, 2) = cos(alpha);
  R(2, 1) = sin(alpha);
  R(1, 2) = -R(2, 1);
  return R;
}

inline Eigen::Matrix3d RZ(const double alpha) {
  Eigen::Matrix3d R = Eigen::Matrix3d::Zero();
  R(2, 2) = 1;
  R(0, 0) = R(1, 1) = cos(alpha);
  R(1, 0) = sin(alpha);
  R(0, 1) = -R(1, 0);
  return R;
}

inline Eigen::Matrix3d RY(const double alpha) {
  Eigen::Matrix3d R = Eigen::Matrix3d::Zero();
  R(1, 1) = 1;
  R(0, 0) = R(2, 2) = cos(alpha);
  R(0, 2) = sin(alpha);
  R(2, 0) = -R(0, 2);
  return R;
}

}

#endif
