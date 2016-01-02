#include "spin_trans.h"

#include <iostream>
#include <jtflib/mesh/util.h>

#include "config.h"
#include "util.h"
#include "arpaca.h"

using namespace std;
using namespace Eigen;
using namespace zjucad::matrix;

namespace riemann {

extern "C"
void conv_quat_to_mat(const double *q, double *m) {
  m[0] = q[0];
  m[1] = q[1];
  m[2] = q[2];
  m[3] = q[3];

  m[4] = -q[1];
  m[5] = q[0];
  m[6] = q[3];
  m[7] = -q[2];

  m[8] = -q[2];
  m[9] = -q[3];
  m[10] = q[0];
  m[11] = q[1];

  m[12] = -q[3];
  m[13] = q[2];
  m[14] = -q[1];
  m[15] = q[0];
}

/// -------------------------------
/// ATTENTION: NODE IS 4xN MATRIX -
/// -------------------------------
spin_trans::spin_trans(const mati_t &tris, const matd_t &nods)
  : tris_(tris), nods_(nods) {
  Mf_.setZero(4*tris_.size(2));
  Mv_.setZero(4*nods_.size(2));
  for (size_t i = 0; i < tris_.size(2); ++i) {
    matd_t vert = nods(colon(1, 3), tris_(colon(), i));
    double area = jtf::mesh::cal_face_area(vert);
    Mf_.segment<4>(4*i) = area*Vector4d::Ones();
    Mv_.segment<4>(4*tris_(0, i)) += 1.0/3*area*Vector4d::Ones();
    Mv_.segment<4>(4*tris_(1, i)) += 1.0/3*area*Vector4d::Ones();
    Mv_.segment<4>(4*tris_(2, i)) += 1.0/3*area*Vector4d::Ones();
  }
}

void spin_trans::set_curvature_change(const matd_t &delta) {
  ASSERT(delta.size() == tris_.size(2));
  rho_ = &delta[0];
}

void spin_trans::build_dirac_operator(SparseMatrix<double> &D) {
  vector<Triplet<double>> trips;
  for (size_t i = 0; i < tris_.size(2); ++i) {
    const double w = -0.5/Mf_[4*i];
    for (size_t j = 0; j < 3; ++j) {
      matd_t ej = nods_(colon(), tris_((j+2)%3, i))-nods_(colon(), tris_((j+1)%3, i));
      matd_t H = zeros<double>(4, 4);
      conv_quat_to_mat(&ej[0], &H[0]);
      H *= w;
      insert_block(4*i, 4*tris_(j, i), &H[0], 4, 4, &trips);
    }
  }
  D.resize(4*tris_.size(2), nods_.size());
  D.reserve(trips.size());
  D.setFromTriplets(trips.begin(), trips.end());
}

void spin_trans::build_rho_operator(SparseMatrix<double> &R) {
  vector<Triplet<double>> trips;
  for (size_t i = 0; i < tris_.size(2); ++i) {
    const double w = rho_[i]/3.0;
    if ( w != 0.0 ) {
      add_diag_block<double, 4>(4*i, 4*tris_(0, i), w, &trips);
      add_diag_block<double, 4>(4*i, 4*tris_(1, i), w, &trips);
      add_diag_block<double, 4>(4*i, 4*tris_(2, i), w, &trips);
    }
  }
  R.resize(4*tris_.size(2), nods_.size());
  R.reserve(trips.size());
  R.setFromTriplets(trips.begin(), trips.end());
}

int spin_trans::solve_eigen_prob() {
  build_dirac_operator(Dirac_);
  build_rho_operator(R_);
  SparseMatrix<double> A = Dirac_-R_;
  SparseMatrix<double> LHS = (A*Mv_.cwiseSqrt().asDiagonal()).transpose()
      *Mf_.asDiagonal()*A*Mv_.cwiseSqrt().asDiagonal();
  if ( is_symm(LHS) ) {
    cerr << "[info] unsymmetric matrix for eigen solve\n";
    return __LINE__;
  }
  arpaca::SymmetricEigenSolver<double> solver
      = arpaca::Solve(LHS, 10, arpaca::MAGNITUDE_SMALLEST);
  lambda_ = solver.eigenvectors().col(0);
  return 0;
}

int spin_trans::solve_poisson_prob(matd_t &nods) {
  Map<VectorXd> X(&nods[0], nods_.size());
  SparseMatrix<double> L;
  VectorXd rhs;
  return 0;
}

}
