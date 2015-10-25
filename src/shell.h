#ifndef DISCRETE_SHELL_H
#define DISCRETE_SHELL_H

#include <zjucad/matrix/matrix.h>

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

namespace riemann {

void get_edge_elem(const mati_t &tris, mati_t &edge);

void get_diamond_elem(const mati_t &tris, mati_t &diam);

class shell_solver
{
public:
  void temp_test() const;
  int solve(double *x) const;
};

}
#endif
