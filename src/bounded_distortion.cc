#include "bounded_distortion.h"

#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>

using namespace std;
using namespace Eigen;
using namespace zjucad::matrix;

namespace riemann {

int calc_tet_base_inv(const mati_t &tets, const matd_t &nods, matd_t &binv) {
  binv.resize(9, tets.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < tets.size(2); ++i) {
    matd_t base = nods(colon(), tets(colon(1, 3), i))-nods(colon(), tets(0, i))*ones<double>(1, 3);
    Map<Matrix3d>(&binv(0, i)) = Map<Matrix3d>(&base[0]).inverse();
  }
  return 0;
}

static inline void AddDiagBlock3d(const size_t i, const size_t j, const double val, vector<Triplet<double>> &trips) {
  trips.push_back(Triplet<double>(3*i+0, 3*j+0, val));
  trips.push_back(Triplet<double>(3*i+1, 3*j+1, val));
  trips.push_back(Triplet<double>(3*i+2, 3*j+2, val));
}

int calc_tet_defo_grad_map(const mati_t &tets, const matd_t &binv, Eigen::SparseMatrix<double> *T) {
  vector<Triplet<double>> trips;
  for (size_t i = 0; i < tets.size(2); ++i) {
    AddDiagBlock3d(3*i+0, tets(0, i), -sum(binv(colon(0, 2), i)), trips);
    AddDiagBlock3d(3*i+0, tets(1, i), binv(0, i), trips);
    AddDiagBlock3d(3*i+0, tets(2, i), binv(1, i), trips);
    AddDiagBlock3d(3*i+0, tets(3, i), binv(2, i), trips);

    AddDiagBlock3d(3*i+1, tets(0, i), -sum(binv(colon(3, 5), i)), trips);
    AddDiagBlock3d(3*i+1, tets(1, i), binv(3, i), trips);
    AddDiagBlock3d(3*i+1, tets(2, i), binv(4, i), trips);
    AddDiagBlock3d(3*i+1, tets(3, i), binv(5, i), trips);

    AddDiagBlock3d(3*i+2, tets(0, i), -sum(binv(colon(6, 8), i)), trips);
    AddDiagBlock3d(3*i+2, tets(1, i), binv(6, i), trips);
    AddDiagBlock3d(3*i+2, tets(2, i), binv(7, i), trips);
    AddDiagBlock3d(3*i+2, tets(3, i), binv(8, i), trips);
  }
  T->resize(9*tets.size(2), 3*(max(tets)+1));
  T->reserve(trips.size());
  T->setFromTriplets(trips.begin(), trips.end());
  return 0;
}

}
