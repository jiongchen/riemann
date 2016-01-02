#ifndef SPIN_TRANS_H
#define SPIN_TRANS_H

#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

namespace riemann {

class spin_trans
{
public:
  spin_trans(const mati_t &tris, const matd_t &nods);
  void set_curvature_change(const matd_t &delta);
  int deform(matd_t &nods);
public:
  void build_dirac_operator(Eigen::SparseMatrix<double> &D);
  void build_rho_operator(Eigen::SparseMatrix<double> &R);
  int solve_eigen_prob();
  int solve_poisson_prob(matd_t &nods);
private:
  const mati_t &tris_;
  const matd_t &nods_;
  const double *rho_;

  Eigen::VectorXd Mf_, Mv_;

  Eigen::VectorXd lambda_;
  Eigen::SparseMatrix<double> Dirac_, R_, L_;
};

}

#endif
