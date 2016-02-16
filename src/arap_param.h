#ifndef ARAP_PARAM_H
#define ARAP_PARAM_H

#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>

namespace riemann {

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;
using cscd_t=Eigen::SparseMatrix<double>;

template <typename T>
class Functional;

class arap_param_solver
{
public:
  arap_param_solver(const mati_t &tris, const matd_t &nods);
  int precompute();
  int solve(double *x0) const;
private:
  std::shared_ptr<Functional<double>> arap_;
  cscd_t LHS_;
  Eigen::SimplicialCholesky<cscd_t> solver_;
};

}

#endif
