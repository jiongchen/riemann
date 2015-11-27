#ifndef BOUNDED_DISTORTION_H
#define BOUNDED_DISTORTION_H

#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>

namespace riemann {

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;
using triplet_t=Eigen::Triplet<double>;

int calc_tet_base_inv(const mati_t &tets, const matd_t &nods, matd_t &binv);
int calc_tet_defo_grad_map(const mati_t &tets, const matd_t &binv, Eigen::SparseMatrix<double> *T);

class bounded_dist_solver
{
public:
  bounded_dist_solver(const mati_t &tets, const matd_t &nods);
private:

};

}

#endif
