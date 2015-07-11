#ifndef GRAD_OPERATOR_H
#define GRAD_OPERATOR_H

#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>

namespace geom_deform {

using mati_t = zjucad::matrix::matrix<size_t>;
using matd_t = zjucad::matrix::matrix<double>;
using spmat_t = Eigen::SparseMatrix<double>;

/**
 * @brief calc_tri_height_vector
 * @param vert: [v0 v1 v2]
 * @param height: [h0 h1 h2]
 *
 * 0------->1
 *          |
 *          |
 *          v
 *          2
 */
template <size_t dim>
void calc_tri_height_vector(const double *vert, double *height) {
  Eigen::Map<const Eigen::Matrix<double, dim, 3>> X(vert);
  Eigen::Map<Eigen::Matrix<double, dim, 3>> H(height);
  Eigen::Matrix<double, dim, dim> I = Eigen::Matrix<double, dim, dim>::Identity();
  for (size_t i = 0; i < 3; ++i) {
    size_t idx[3] = {i, (i+1)%3, (i+2)%3};
    Eigen::Matrix<double, dim, 1> x12 = X.col(idx[2])-X.col(idx[1]);
    H.col(i) = (I-x12*x12.transpose()/x12.dot(x12))*(X.col(idx[0])-X.col(idx[1]));
  }
}

/**
 * @brief calc_grad_operator
 * @param tris
 * @param nods
 * @param G: #face*3 by #vert, G*u, u is a scalar field
 *        defined on vertex
 */
void calc_grad_operator(const mati_t &tris, const matd_t &nods, spmat_t *G);

}

#endif
