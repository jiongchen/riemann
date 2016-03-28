#include "conformal_volume.h"

#include <iostream>
#include <fstream>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/lapack.h>

#include "config.h"
#include "grad_operator.h"
#include "util.h"
#include "vtk.h"

using namespace std;
using namespace Eigen;
using namespace zjucad::matrix;

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

template <size_t dim>
static void mass_matrix(const mati_t &tets, const matd_t &nods, SparseMatrix<double> *M) {
  vector<Triplet<double>> trips;
  for (size_t i = 0; i < tets.size(2); ++i) {
    matd_t ed = nods(colon(), tets(colon(1, 3), i))-nods(colon(), tets(0, i))*ones<double>(1, 3);
    double vol = fabs(det(ed))/6.0;
    add_diag_block<double, dim>(tets(0, i), tets(0, i), vol/4.0, &trips);
    add_diag_block<double, dim>(tets(1, i), tets(1, i), vol/4.0, &trips);
    add_diag_block<double, dim>(tets(2, i), tets(2, i), vol/4.0, &trips);
    add_diag_block<double, dim>(tets(3, i), tets(3, i), vol/4.0, &trips);
  }
  const size_t n = nods.size(2);
  M->resize(dim*n, dim*n);
  M->setFromTriplets(trips.begin(), trips.end());
}

template <size_t dim>
static void laplacian_matrix(const mati_t &tets, const matd_t &nods, SparseMatrix<double> *L) {
  vector<Triplet<double>> trips;
  for (size_t i = 0; i < tets.size(2); ++i) {
    matd_t v = nods(colon(), tets(colon(), i));
    matd_t g = zeros<double>(3, 4);
    calc_tet_linear_basis_grad(&v[0], &g[0]);
    matd_t ed = v(colon(), colon(1, 3))-v(colon(), 0)*ones<double>(1, 3);
    double vol = fabs(det(ed))/6.0;
    for (size_t p = 0; p < 4; ++p) {
      for (size_t q = 0; q < 4; ++q) {
        double val = dot(g(colon(), p), g(colon(), q))*vol;
        add_diag_block<double, dim>(tets(p, i), tets(q, i), val, &trips);
      }
    }
  }
  const size_t n = nods.size(2);
  L->resize(dim*n, dim*n);
  L->setFromTriplets(trips.begin(), trips.end());
}

//==============================================================================
conformal_volume::conformal_volume(const mati_t &tets, const matd_t &verts)
  : tets_(tets), verts_(verts) {
  u_.setZero(verts_.size(2));
  gradu_.resize(NoChange, tets_.size(2));
  gradu_.setZero();
  vol_ = zeros<double>(tets_.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < tets_.size(2); ++i) {
    matd_t e = verts_(colon(1, 3), tets_(colon(1, 3), i))-verts_(colon(1, 3), tets_(0, i))*ones<double>(1, 3);
    vol_[i] = fabs(det(e))/6.0;
  }
  I_ << 0, -1, 0, 0,
      1, 0, 0, 0,
      0, 0, 0, -1,
      0, 0, 1, 0;
  J_ << 0, 0, -1, 0,
      0, 0, 0, 1,
      1, 0, 0, 0,
      0, -1, 0, 0;
  K_ << 0, 0, 0, -1,
      0, 0, -1, 0,
      0, 1, 0, 0,
      1, 0, 0, 0;
}

void conformal_volume::set_charge(const double *pos, const double intensity) {
  ASSERT(intensity > 0.0);
  Vector3d a(pos);
  // ->FOR SCALAR
  VectorXd psi(verts_.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < psi.size(); ++i) {
    Vector3d f(verts_(1, i), verts_(2, i), verts_(3, i));
    psi(i) = intensity/(f-a).norm();
  }
//  psi.cwise().log()/2;
  // ->FOR GRADIENT
  Matrix3Xd gu;
  calc_grad_u(psi, gu);
  // ->APPEND
  u_ += psi;
  gradu_.bottomLeftCorner(3, gradu_.cols()) += gu;
}

void conformal_volume::calc_grad_u(const VectorXd &u, Matrix3Xd &grad) {
  grad.resize(NoChange, tets_.size(2));
  itr_matrix<const double *> U(u.size(), 1, u.data());
  itr_matrix<double *> G(grad.rows(), grad.cols(), grad.data());
#pragma omp parallel for
  for (size_t i = 0; i < tets_.size(2); ++i) {
    matd_t v = verts_(colon(1, 3), tets_(colon(), i));
    matd_t g = zeros<double>(3, 4);
    calc_tet_linear_basis_grad(&v[0], &g[0]);
    G(colon(), i) = g*U(tets_(colon(), i));
  }
}

void conformal_volume::solve_eigen_prob() {
  SparseMatrix<double> E;
  laplacian_matrix<4>(tets_, verts_(colon(1, 3), colon()), &L_);
  // get sqrM
  VectorXd sqrinvM(4*verts_.size(2)); {
    for (size_t i = 0; i < tets_.size(2); ++i)
      for (size_t j = 0; j < 4; ++j)
        sqrinvM.segment<4>(4*tets_(j, i)) += vol_[i]/4.0*Vector4d::Ones();
    sqrinvM.cwiseSqrt().cwiseInverse();
  }
  // assemble B
  SparseMatrix<double> B; {
    vector<Triplet<double>> trips;
    for (size_t i = 0; i < tets_.size(2); ++i) {
      matd_t v = verts_(colon(1, 3), tets_(colon(), i));
      Matrix<double, 3, 4> g = Matrix<double, 3, 4>::Zero();
      calc_tet_linear_basis_grad(&v[0], g.data());
      for (size_t p = 0; p < 4; ++p) {
        for (size_t q = 0; q < 4; ++q) {
          Vector4d XG = 1.0/12*(I_*gradu_.col(i)*3*vol_[i]*g.col(q).dot(Vector3d(1, 0, 0))
                                +J_*gradu_.col(i)*3*vol_[i]*g.col(q).dot(Vector3d(0, 1, 0))
                                +K_*gradu_.col(i)*3*vol_[i]*g.col(q).dot(Vector3d(0, 0, 1)));
          Matrix4d mat;
          conv_quat_to_mat(XG.data(), mat.data());
          insert_block(4*tets_(p, i), 4*tets_(q, i), mat.data(), 4, 4, &trips);
        }
      }
    }
    B.resize(4*verts_.size(2), 4*verts_.size(2));
    B.setFromTriplets(trips.begin(), trips.end());
  }
  // assemble weighted M
  SparseMatrix<double> gM; {
    vector<Triplet<double>> trips;
    for (size_t i = 0; i < tets_.size(2); ++i) {
      double w = 0.75*gradu_.col(i).squaredNorm();
      for (size_t p = 0; p < 4; ++p) {
        for (size_t q = 0; q < 4; ++q) {
          double m = (p == q) ? vol_[i]/10 : vol_[i]/20;
          add_diag_block(tets_(p, i), tets_(q, i), w*m, &trips);
        }
      }
    }
    gM.resize(4*verts_.size(2), 4*verts_.size(2));
    gM.setFromTriplets(trips.begin(), trips.end());
  }
  SparseMatrix<double> BT = B.transpose();
  E = L_+0.5*(B+BT)+gM;
  SparseMatrix<double> LHS = sqrinvM.asDiagonal()*E*sqrinvM.asDiagonal();

  // inverse power iteration
  SimplicialCholesky<SparseMatrix<double>> solver;
  solver.compute(LHS);
  ASSERT(solver.info() == Success);
  VectorXd eigvec = VectorXd::Ones(LHS.cols());
  eigvec /= eigvec.norm();
  for (size_t iter = 0; iter < 3; ++iter) {
    eigvec = (solver.solve(eigvec)).eval();
    ASSERT(solver.info() == Success);
    eigvec /= eigvec.norm();
  }
  lambda_ = sqrinvM.asDiagonal()*eigvec;
}

void conformal_volume::solve_poisson_prob(double *x) {
  Map<VectorXd> X(x, 4*verts_.size(2));
  VectorXd rhs; { // calculate the divergence

  }
  SimplicialCholesky<SparseMatrix<double>> solver;
  solver.compute(L_);
  ASSERT(solver.info() == Success);
  X = solver.solve(rhs);
  ASSERT(solver.info() == Success);
}

//==============================================================================
void conformal_volume::draw_gradient(const char *filename) {
  mati_t line(2, tets_.size(2));
  line(0, colon()) = colon(0, tets_.size(2)-1);
  line(1, colon()) = colon(tets_.size(2), 2*tets_.size(2)-1);
  matd_t nods(3, 2*tets_.size(2));
  itr_matrix<const double *> G(gradu_.rows(), gradu_.cols(), gradu_.data());
  for (size_t i = 0; i < tets_.size(2); ++i) {
    nods(colon(), i) = verts_(colon(1, 3), tets_(colon(), i))*ones<double>(4, 1)/4;
    nods(colon(), i+tets_.size(2)) = nods(colon(), i)+0.02*G(colon(1, 3), i);
  }
  ofstream os(filename);
  line2vtk(os, &nods[0], nods.size(2), &line[0], line.size(2));
  os.close();
}

void conformal_volume::debug_laplacian() {
  u_.setZero();
  u_(6) = 100.0;
  u_(608) = -100.0;
  vector<size_t> g2l(u_.size());
  size_t cnt = 0;
  for (size_t i = 0; i < g2l.size(); ++i) {
    if ( i == 6 || i == 608 )
      g2l[i] = -1;
    else
      g2l[i] = cnt++;
  }
  laplacian_matrix<1>(tets_, verts_(colon(1, 3), colon()), &L_);
  VectorXd rhs = VectorXd::Zero(u_.size())-L_*u_;
  rm_spmat_col_row(L_, g2l);
  rm_vector_row(rhs, g2l);
  SimplicialCholesky<SparseMatrix<double>> solver;
  solver.compute(L_);
  ASSERT(solver.info() == Eigen::Success);
  VectorXd du = solver.solve(rhs);
  ASSERT(solver.info() == Eigen::Success);
  VectorXd DU(u_.size());
  rc_vector_row(du, g2l, DU);
  u_ += DU;
}

}
