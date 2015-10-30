#ifndef DISCRETE_SHELL_H
#define DISCRETE_SHELL_H

#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

namespace riemann {

template<typename T>
class Functional;

template<typename T>
class Constraint;

struct shell_args {
  double ws, wb, wp;
  size_t max_iter;
  double tolerance;
};

class shell_deformer
{
public:
  shell_deformer(const mati_t &tris, const matd_t &nods, shell_args &args);
  void fix_vert(const size_t id, const double *coords);
  void free_vert(const size_t id);
  int prepare();
  int solve(double *x);
  void unit_test() const;
private:
  shell_args args_;
  mati_t edges_, diams_;
  std::vector<std::shared_ptr<Constraint<double>>> cbf_;
  std::shared_ptr<Constraint<double>> constraint_;
  Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver_;
};

}
#endif
