#ifndef CONFORMAL_VOLUME_H
#define CONFORMAL_VOLUME_H

#include <Eigen/Sparse>
#include <zjucad/matrix/matrix.h>

namespace riemann {

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

class conformal_volume
{
public:
  conformal_volume(const mati_t &tets, const matd_t &verts);
  void set_charge(const double *pos, const double intensity);
  void solve_eigen_prob();
  void solve_poisson_prob(double *x);
  // DEBUG
  void draw_gradient(const char *filename);
  void debug_laplacian();
private:
  void calc_grad_u(const Eigen::VectorXd &u, Eigen::Matrix3Xd &grad);
public:
  const mati_t &tets_;
  const matd_t &verts_;

  Eigen::VectorXd lambda_;
  Eigen::VectorXd u_;
  Eigen::Matrix4Xd gradu_;

  Eigen::SparseMatrix<double> L_, M_, B_;
};

}

#endif
