#ifndef HOLE_FILLING_H
#define HOLE_FILLING_H

#include <Eigen/Dense>

namespace riemann {

using mati=Eigen::Matrix<int, -1, -1, Eigen::RowMajor>;
using matd=Eigen::Matrix<double, -1, -1, Eigen::RowMajor>;
using veci=Eigen::Matrix<int, -1, 1>;
using vecd=Eigen::Matrix<double, -1, 1>;

struct boundary_loop {
  veci bnd;
  matd pos;
  vecd sf;
  matd T, B, N;
  vecd dl;
};

void get_boundary_loop(const mati &tris, const matd &nods, boundary_loop &loop);
void calc_bnd_local_frame(const mati &tris, const matd &nods, boundary_loop &loop);
void calc_dl(boundary_loop &loop);
void init_bnd_scalar_field(boundary_loop &loop);

double indicator_value(const vecd &x, const boundary_loop &loop);
void uniform_sampling(const matd &box, const int num, matd &pts);
//void calc_scalar_field(const boundary_loop &loop, const MatrixXd &pts, Eigen::VectorXd &sf);

}

#endif
