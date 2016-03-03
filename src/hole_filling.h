#ifndef HOLE_FILLING_H
#define HOLE_FILLING_H

#include <Eigen/Dense>

namespace riemann {

using eigen_mati_t=Eigen::Matrix<int, -1, -1, Eigen::RowMajor>;
using eigen_matd_t=Eigen::Matrix<double, -1, -1, Eigen::RowMajor>;
using eigen_veci_t=Eigen::Matrix<int, 1, -1, Eigen::RowMajor>;
using eigen_vecd_t=Eigen::Matrix<double, 1, -1, Eigen::RowMajor>;

struct boundary_loop {
  enum type {OFFSET, SCALE};
  eigen_veci_t bnd;
  eigen_matd_t pos;
  eigen_vecd_t sf;
  eigen_matd_t T, B, N;
  eigen_vecd_t dl;
  int method;
};

void get_boundary_loop(const eigen_mati_t &tris, const eigen_matd_t &nods, boundary_loop &loop);
void select_method(boundary_loop &loop, boundary_loop::type t);
void calc_bnd_local_frame(const eigen_mati_t &tris, const eigen_matd_t &nods, boundary_loop &loop);
void calc_infinitesimal_elem(boundary_loop &loop);
void calc_bnd_indicator(boundary_loop &loop);
void uniform_random_sampling(const eigen_matd_t &box, const int num, eigen_matd_t &pts);
void structured_grid_sampling(const eigen_matd_t &box, const int res, eigen_matd_t &pts);
double indicator_value(const eigen_vecd_t &x, const boundary_loop &loop);
void calc_scalar_field(const boundary_loop &loop, const eigen_matd_t &pts, eigen_vecd_t &sf);

}

#endif
