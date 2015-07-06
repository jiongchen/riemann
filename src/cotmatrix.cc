#include "cotmatrix.h"

#include <set>

#include "util.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

namespace surfparam {

/**
 * Discrete Laplace-Beltrami Operator for Triangle Mesh
 * (\triangle f)(v_i) = \sum_{j \in \mathcal{N}_i} w_{ij}(f(v_j)-f(v_i))
 */

template <typename T>
void compute_vert_valence(const matrix<size_t> &cell,
                          const matrix<double> &nods,
                          const size_t dim,
                          Matrix<T, -1, 1> &valence) {
  set<pair<size_t, size_t>> vis;
  valence.setZero(dim*nods.size(2), 1);
  Matrix<T, -1, 1> ones;
  ones.setOnes(dim, 1);
  for (size_t i = 0; i < cell.size(2); ++i) {
    for (size_t j = 0; j < 3; ++j) {
      size_t vp = cell(j, i);
      size_t vq = cell((j+1)%3, i);
      auto edge = (vp < vq ? make_pair(vp, vq) : make_pair(vq, vp));
      if ( vis.find(edge) != vis.end() )
        continue;
      vis.insert(edge);
      valence.segment(dim*vp, dim) += ones;
      valence.segment(dim*vq, dim) += ones;
    }
  }
}

/**
 * @brief w_{ij} = \frac{1}{2}(cot \alpha_{ij} + cot \beta_{ij})
 */
void cotmatrix(const matrix<size_t> &cell,
               const matrix<double> &nods,
               const size_t dim,
               SparseMatrix<double> *L,
               bool normalized) {
  const size_t edge[3][2] = {{1, 2}, {2, 0}, {0, 1}};
  std::vector<Triplet<double>> trips;

  for (size_t i = 0; i < cell.size(2); ++i) {
    matrix<double> vert = nods(colon(), cell(colon(), i));
    matrix<double> half_cot_val(3);
    half_cot_val[0] = 0.5*cal_cot_val(&vert(0, 1), &vert(0, 0), &vert(0, 2));
    half_cot_val[1] = 0.5*cal_cot_val(&vert(0, 0), &vert(0, 1), &vert(0, 2));
    half_cot_val[2] = 0.5*cal_cot_val(&vert(0, 0), &vert(0, 2), &vert(0, 1));

    for (size_t k = 0; k < 3; ++k) {
      size_t src = cell(edge[k][0], i);
      size_t des = cell(edge[k][1], i);
      runtime_dim_add_diag_block(dim, src, src, -half_cot_val[k], &trips);
      runtime_dim_add_diag_block(dim, des, des, -half_cot_val[k], &trips);
      runtime_dim_add_diag_block(dim, src, des, half_cot_val[k], &trips);
      runtime_dim_add_diag_block(dim, des, src, half_cot_val[k], &trips);
    }
  }
  const size_t lap_size = dim * nods.size(2);
  L->resize(lap_size, lap_size);
  L->reserve(trips.size());
  L->setFromTriplets(trips.begin(), trips.end());

  if ( normalized ) {

  }
}

/**
 * @brief w_{ij} = 1
 */
void unimatrix(const matrix<size_t> &cell,
               const matrix<double> &nods,
               const size_t dim,
               SparseMatrix<double> *L,
               bool normalized) {
  vector<Triplet<double>> trips;
  set<pair<size_t, size_t>> vis;

  for (size_t i = 0; i < cell.size(2); ++i) {
    for (size_t j = 0; j < 3; ++j) {
      size_t src = cell(j, i);
      size_t des = cell((j+1)%3, i);
      auto edge = (src < des ? make_pair(src, des) : make_pair(des, src));
      if ( vis.find(edge) != vis.end() )
        continue;
      vis.insert(edge);
      runtime_dim_add_diag_block(dim, src, src, -1.0, &trips);
      runtime_dim_add_diag_block(dim, des, des, -1.0, &trips);
      runtime_dim_add_diag_block(dim, src, des, 1.0, &trips);
      runtime_dim_add_diag_block(dim, des, src, 1.0, &trips);
    }
  }
  const size_t lap_size = dim * nods.size(2);
  L->resize(lap_size, lap_size);
  L->reserve(trips.size());
  L->setFromTriplets(trips.begin(), trips.end());

  if ( normalized ) {
    VectorXd valence;
    compute_vert_valence(cell, nods, dim, valence);
    *L = (valence.asDiagonal().inverse()*(*L)).eval();
  }
}

}
