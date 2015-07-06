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
 */
void calc_tri_height_vector(const double *vert, double *height);

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
