#ifndef PARAM_UTIL_H
#define PARAM_UTIL_H

#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <unordered_set>

#include "config.h"

namespace surfparam {

template <typename T, size_t dim = 3>
void add_diag_block(const size_t row, const size_t col, const T val, std::vector<Eigen::Triplet<T>> *mat) {
  const size_t row_start = row * dim;
  const size_t col_start = col * dim;
  for (size_t offset = 0; offset < dim; ++offset)
    mat->push_back(Eigen::Triplet<T>(row_start+offset, col_start+offset, val));
}

template <typename T>
void runtime_dim_add_diag_block(const size_t dim, const size_t row, const size_t col, const T val,
                                std::vector<Eigen::Triplet<T>> *trip) {
  const size_t row_start = row * dim;
  const size_t col_start = col * dim;
  for (size_t off = 0; off < dim; ++off)
    trip->push_back(Eigen::Triplet<T>(row_start+off, col_start+off, val));
}

template <typename T>
void rm_spmat_col_row(Eigen::SparseMatrix<T> &A,
                      const std::unordered_set<size_t> &idx) {
  std::vector<size_t> g2l(A.cols());
  size_t ptr = 0;
  for (size_t i = 0; i < g2l.size(); ++i) {
    if ( idx.find(i) != idx.end() )
      g2l[i] = -1;
    else
      g2l[i] = ptr++;
  }
  rm_spmat_col_row<T>(A, g2l);
}

template <typename T>
void rm_spmat_col_row(Eigen::SparseMatrix<T> &A,
                      const std::vector<size_t> &g2l) {
  size_t new_size = 0;
  for (size_t i = 0; i < g2l.size(); ++i) {
    if ( g2l[i] != -1)
      ++new_size;
  }
  std::vector<Eigen::Triplet<T>> trips;
  for (size_t j = 0; j < A.outerSize(); ++j) {
    for (typename Eigen::SparseMatrix<T>::InnerIterator it(A, j); it; ++it) {
      if ( g2l[it.row()] != -1 && g2l[it.col()] != -1 )
        trips.push_back(Eigen::Triplet<T>(g2l[it.row()], g2l[it.col()], it.value()));
    }
  }
  A.resize(new_size, new_size);
  A.reserve(trips.size());
  A.setFromTriplets(trips.begin(), trips.end());
}

template <typename T>
void rm_vector_row(Eigen::Matrix<T, -1, 1> &b,
                   const std::vector<size_t> &g2l) {
  size_t new_size = 0;
  for (size_t i = 0; i < g2l.size(); ++i) {
    if ( g2l[i] != -1 )
      ++new_size;
  }
  Eigen::Matrix<T, -1, 1> sub;
  sub.resize(new_size);
#pragma omp parallel for
  for (size_t i = 0; i < g2l.size(); ++i)
    if ( g2l[i] != -1 )
      sub[g2l[i]] = b[i];
  b = sub;
}

template <typename T>
void up_vector_row(const Eigen::Matrix<T, -1, 1> &l, const std::vector<size_t> &g2l, Eigen::Matrix<T, -1, 1> &g) {
#pragma omp parallel for
  for (size_t i = 0; i < g2l.size(); ++i) {
    if ( g2l[i] != -1 )
      g[i] = l[g2l[i]];
  }
}

template <typename T>
Eigen::SparseMatrix<T> sparse_diag_matrix(const Eigen::DiagonalMatrix<T, -1> &diag) {
  std::vector<Eigen::Triplet<T>> trips;
  for (size_t i = 0; i < diag.cols(); ++i)
    trips.push_back(Eigen::Triplet<T>(i, i, diag.diagonal()[i]));
  Eigen::SparseMatrix<T> rtn;
  rtn.resize(diag.rows(), diag.cols());
  rtn.reserve(trips.size());
  rtn.setFromTriplets(trips.begin(), trips.end());
  return rtn;
}

template <typename T>
bool is_symm(const Eigen::SparseMatrix<T> &A) {
  Eigen::SparseMatrix<T> AT = A.transpose();
  if ( (AT - A).norm() < 1e-16 )
    return true;
  return false;
}

template <typename T>
inline T cal_cot_val(const T* a, const T* b, const T* c) {
  Eigen::Matrix<T, 3, 1> ab(b[0]-a[0], b[1]-a[1], b[2]-a[2]);
  Eigen::Matrix<T, 3, 1> bc(c[0]-b[0], c[1]-b[1], c[2]-b[2]);
  Eigen::Matrix<T, 3, 1> ca(a[0]-c[0], a[1]-c[1], a[2]-c[2]);
  return 0.5 * (ab.dot(ab) + bc.dot(bc) - ca.dot(ca)) / ab.cross(bc).norm();
}

template <typename T>
inline T cal_cot_val(const T *point_set) {
  return cal_cot_val<T>(&point_set[0], &point_set[3], &point_set[6]);
}

}
#endif
