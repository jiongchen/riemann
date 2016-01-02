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
  void precompute();
  int deform(matd_t &x);
private:
  void build_dirac_operator(Eigen::SparseMatrix<double> &D);
  void build_rho_operator(Eigen::SparseMatrix<double> &R);
  void calc_div_f(const matd_t &x, Eigen::VectorXd divf);
  int solve_eigen_prob(Eigen::VectorXd &lambda);
  int solve_poisson_prob(const Eigen::VectorXd &lambda, matd_t &x);
private:
  const mati_t &tris_;
  const matd_t &nods_;
  matd_t rho_;

  Eigen::VectorXd Mf_, Mv_;
  Eigen::SparseMatrix<double> Dirac_, R_, L_;
  Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> ldlt_solver_;
};

}

#endif
