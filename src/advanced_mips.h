#ifndef ADVANCED_MIPS_H
#define ADVANCED_MIPS_H

#include <unordered_set>
#include <zjucad/matrix/matrix.h>

namespace riemann {

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

class move_vertex;

class mips_deformer_2d
{
public:
  mips_deformer_2d(const mati_t &tris, const matd_t &nods);
  int deform(double *x, const size_t maxiter=20000) const;
  void unit_test() const;
private:
  int apply(double *x) const;
private:
  const size_t dim_;
  std::unordered_set<size_t> fixed_;
  std::shared_ptr<move_vertex> move_;
};

}

#endif
