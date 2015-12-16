#ifndef ADVANCED_MIPS_H
#define ADVANCED_MIPS_H

#include <zjucad/matrix/matrix.h>

namespace riemann {

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

class mips_deformer_2d
{
public:
  mips_deformer_2d(const mati_t &tris, const matd_t &nods);
  int deform(double *x) const;
  void unit_test() const;
private:

};

}

#endif
