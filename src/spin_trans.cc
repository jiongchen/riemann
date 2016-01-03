#include "spin_trans.h"

#include <iostream>
#include <jtflib/mesh/util.h>
#include <jtflib/mesh/mesh.h>
#include <zjucad/matrix/io.h>

#include "config.h"
#include "util.h"
#include "arpaca.h"
#include "cotmatrix.h"

using namespace std;
using namespace Eigen;
using namespace zjucad::matrix;
using namespace jtf::mesh;

namespace riemann {

static inline void conv_quat_to_mat(const double *q, double *m) {
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

static inline Vector4d conjugate(const Vector4d &q) {
  return Vector4d(q[0], -q[1], -q[2], -q[3]);
}

static inline Vector4d quat_prod(const Vector4d &a, const Vector4d &b) {
  double w = a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3],
      x = a[2]*b[3]-a[3]*b[2]+a[0]*b[1]+a[1]*b[0],
      y = -a[1]*b[3]+a[0]*b[2]+a[3]*b[1]+a[2]*b[0],
      z = a[0]*b[3]+a[1]*b[2]-a[2]*b[1]+a[3]*b[0];
  return Vector4d(w, x, y, z);
}

/// -------------------------------
/// ATTENTION: NODE IS 4xN MATRIX -
/// -------------------------------
spin_trans::spin_trans(const mati_t &tris, const matd_t &nods)
  : tris_(tris), nods_(nods) {
  ASSERT(nods_.size(1) == 4);
  Mf_.setZero(4*tris_.size(2));
  Mv_.setZero(4*nods_.size(2));
  for (size_t i = 0; i < tris_.size(2); ++i) {
    matd_t vert = nods_(colon(1, 3), tris_(colon(), i));
    const double area = jtf::mesh::cal_face_area(vert);
    Mf_.segment<4>(4*i) = area*Vector4d::Ones();
    Mv_.segment<4>(4*tris_(0, i)) += area/3*Vector4d::Ones();
    Mv_.segment<4>(4*tris_(1, i)) += area/3*Vector4d::Ones();
    Mv_.segment<4>(4*tris_(2, i)) += area/3*Vector4d::Ones();
  }
  rho_ = zeros<double>(tris_.size(2), 1);
  cotmatrix(tris_, nods_(colon(1, 3), colon()), 4, &L_);
  L_ *= -1; // make it semi-positive definite
  build_dirac_operator(Dirac_);
  ASSERT(Dirac_.cols() == L_.cols());
}

void spin_trans::set_curvature_change(const matd_t &delta) {
  ASSERT(delta.size() == tris_.size(2));
  rho_ = delta;
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
      add_diag_block<double, 4>(i, tris_(0, i), w, &trips);
      add_diag_block<double, 4>(i, tris_(1, i), w, &trips);
      add_diag_block<double, 4>(i, tris_(2, i), w, &trips);
    }
  }
  R.resize(4*tris_.size(2), nods_.size());
  R.reserve(trips.size());
  R.setFromTriplets(trips.begin(), trips.end());
}

int spin_trans::solve_eigen_prob(VectorXd &lambda) {
  /// solve integrable condition: $(D-\rho)\lambda=0$, given
  /// arbitrary $\rho$, we solve an eigenvalue problem to find
  /// the smallest shift s.t. $(D-\rho)\lambda=\gamma\lambda$
  cout << "[INFO] solve for similarity field\n";
  build_rho_operator(R_);
  ASSERT(R_.cols() == L_.cols());
  SparseMatrix<double> A = Dirac_-R_;
  VectorXd V = Mv_.cwiseSqrt().cwiseInverse();
  SparseMatrix<double> LHS = (A*V.asDiagonal()).transpose()*Mf_.asDiagonal()*A*V.asDiagonal();
  if ( !is_symm(LHS) ) {
    cerr << "[INFO] unsymmetric matrix for eigen solve\n";
    return __LINE__;
  }
  arpaca::SymmetricEigenSolver<double> solver;
  solver.SetEigenvalueType(arpaca::MAGNITUDE_SMALLEST);
  solver.SetMaxIterations(20000);
  solver.SetNumLanczosVectors(0);
  solver.SetTolerance(-1);
  solver.Solve(LHS.cols(), 10, arpaca::MakeDefaultOperator(LHS));
  printf("[INFO] arpack %d iter, %d converged, %s\n",
         solver.num_actual_iterations(), solver.num_converged_eigenvalues(), solver.GetInfo());
  ASSERT(V.size() == solver.eigenvectors().rows());
  cout << "[INFO] global constant shift: " << solver.eigenvalues()[0] << endl;
  lambda = V.asDiagonal()*solver.eigenvectors().col(0);
  return 0;
}

void spin_trans::calc_div_f(const VectorXd &lambda, VectorXd &divf) {
  divf.setZero(nods_.size());
  for (size_t i = 0; i < tris_.size(2); ++i) {
    matd_t vert = nods_(colon(), tris_(colon(), i));
    for (size_t j = 0; j < 3; ++j) {
      const size_t p = j, q = (j+1)%3;
      matd_t e = vert(colon(), q)-vert(colon(), p);
      Vector4d epq(&e[0]), lp(&lambda[4*tris_(p, i)]), lq(&lambda[4*tris_(q, i)]),
          lpc = conjugate(lp), lqc = conjugate(lq);
      Vector4d etilde =
           1.0/3*quat_prod(quat_prod(lpc, epq), lp)
          +1.0/6*quat_prod(quat_prod(lpc, epq), lq)
          +1.0/6*quat_prod(quat_prod(lqc, epq), lp)
          +1.0/3*quat_prod(quat_prod(lqc, epq), lq);
      double cota = cal_cot_val(&vert(1, p), &vert(1, (j+2)%3), &vert(1, q));
      divf.segment<4>(4*tris_(p, i)) -= 0.5*cota*etilde;
      divf.segment<4>(4*tris_(q, i)) += 0.5*cota*etilde;
    }
  }
}

int spin_trans::solve_poisson_prob(const VectorXd &lambda, matd_t &x) {
  ASSERT(x.size(1) == nods_.size(1) && x.size(2) == nods_.size(2));
  Map<VectorXd> X(&x[0], x.size());
  VectorXd divf;
  calc_div_f(lambda, divf);
  divf *= -1.0; // corresponds to the change on Laplacian

  // fix the first point to remove the kernel of Laplacian
  vector<size_t> g2l(L_.cols());
  size_t cnt = 0;
  for (size_t i = 0; i < g2l.size(); ++i) {
    if ( i/4 == 0 )
      g2l[i] = -1;
    else
      g2l[i] = cnt++;
  }
  VectorXd rhs = divf-L_*X;
  rm_spmat_col_row(L_, g2l);
  rm_vector_row(rhs, g2l);

  solver_.setMode(SimplicialCholeskyLLT);
  solver_.compute(L_);
  ASSERT(solver_.info() == Success);
  VectorXd dx = solver_.solve(rhs);
  ASSERT(solver_.info() == Success);
  VectorXd Dx = VectorXd::Zero(x.size());
  rc_vector_row(dx, g2l, Dx);
  X += Dx;
  return 0;
}

int spin_trans::deform(matd_t &x) {
  VectorXd lambda;
  solve_eigen_prob(lambda);
  solve_poisson_prob(lambda, x);
  return 0;
}

}
