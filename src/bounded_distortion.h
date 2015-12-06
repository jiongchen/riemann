#ifndef BOUNDED_DISTORTION_H
#define BOUNDED_DISTORTION_H

#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>

namespace riemann {

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;
using triplet_t=Eigen::Triplet<double>;

int calc_tet_base_inv(const mati_t &tets, const matd_t &nods, matd_t &binv);
int calc_tet_df_map(const mati_t &tets, const matd_t &binv, Eigen::SparseMatrix<double> *T);

class bd_pos_constraint;

struct bd_args {
  size_t maxiter;
  double tolerance;
};

class bd_solver
{
public:
  bd_solver(const mati_t &tets, const matd_t &nods, const bd_args &args);
  void set_bound(const double K);
  int pin_down_vert(const size_t id, const double *pos);
  int prefactorize();
  int solve(double *initX) const;
  int alter_solve(double *initX) const;
private:
  int euclidean_proj(const double *Tx, double *PTx) const;

  const mati_t &tets_;
  const matd_t &nods_;
  const size_t dim_, lift_dim_;
  double K_;
  bd_args args_;

  Eigen::SparseMatrix<double> T_;
  Eigen::VectorXd volume_;
  std::shared_ptr<bd_pos_constraint> linc_;
};

}

#endif
