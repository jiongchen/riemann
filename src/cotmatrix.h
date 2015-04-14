#ifndef COTMATRIX_H
#define COTMATRIX_H

#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>

namespace surfparam {

void cotmatrix(const zjucad::matrix::matrix<size_t> &cell,
               const zjucad::matrix::matrix<double> &nods,
               const size_t dim,
               Eigen::SparseMatrix<double> *L);

}
#endif
