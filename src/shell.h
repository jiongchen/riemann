#ifndef DISCRETE_SHELL_H
#define DISCRETE_SHELL_H

#include <zjucad/matrix/matrix.h>

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

namespace riemann {

template<typename T>
class Functional;

template<typename T>
class Constraint;

class shell_deformer
{
public:
  shell_deformer(const mati_t &tris, const matd_t &nods);
  void temp_test() const;
  int solve(double *x) const;
private:
  const mati_t &tris_;
  const matd_t &nods_;
  std::vector<std::shared_ptr<Functional<double>>> ebf_;
  std::shared_ptr<Functional<double>> energy_;
  std::vector<std::shared_ptr<Constraint<double>>> cbf_;
  std::shared_ptr<Constraint<double>> constraint_;
};

}
#endif
