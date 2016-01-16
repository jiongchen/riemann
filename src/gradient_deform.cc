#include "gradient_deform.h"

#include <zjucad/matrix/itr_matrix.h>
#include <jtflib/mesh/io.h>
#include <Eigen/Geometry>

#include "cotmatrix.h"
#include "grad_operator.h"
#include "vtk.h"
#include "util.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

namespace riemann {

static void calc_tris_area(const mati_t &tris, const matd_t &nods, VectorXd &area) {
  area.setZero(tris.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < tris.size(2); ++i) {
    matd_t edge = nods(colon(), tris(colon(1, 2), i))-nods(colon(), tris(colon(0, 1), i));
    area[i] = norm(cross(edge(colon(), 0), edge(colon(), 1)))/2.0;
  }
}

gradient_field_deform::gradient_field_deform(const mati_t &tris, const matd_t &nods)
  : tris_(tris), nods_(nods) {
  // laplacian
  cotmatrix(tris, nods_, 1, &L_);
  L_ *= -1.0;
  // discrete gradient operator
  calc_grad_operator(tris_, nods_, &G_);
  // gradient
  Gxyz_.setZero(3*tris_.size(2), nods_.size(2)); {
    matd_t nodsT = trans(nods);
    Gxyz_ = G_*Map<const MatrixXd>(&nodsT[0], nodsT.size(1), nodsT.size(2));
  }
  // area
  calc_tris_area(tris_, nods_, area_);
  // harmonic field
  hf_.setZero(nods_.size(2));
  // transform field
  // ...
}

//int gradient_field_deform::calc_divergence(const VectorXd &vf, VectorXd &div) const {
//  div.setZero(nods_.size(2));
//  for (size_t i = 0; i < tris_.size(2); ++i) {
//    for (size_t j = 0; j < 3; ++j) {
//      div[tris_(j, i)] += gradB_.col(3*i+j).dot(vf.segment<3>(3*i))*area_[i];
//    }
//  }
//  return 0;
//}

void gradient_field_deform::set_fixed_verts(const vector<size_t> &idx) {
  fixDoF_.clear();
  for (auto &id : idx)
    fixDoF_.insert(id);
}

void gradient_field_deform::set_edited_verts(const vector<size_t> &idx) {
  editDoF_.clear();
  for (auto &id : idx)
    editDoF_.insert(id);
}

//int gradient_field_deform::edit_boundary(const vector<size_t> &idx) {
//  editDOF_.clear();
//  for (auto &id : idx) {
//    editDOF_.insert(id);
//    hf_[id] = 1.0;
//  }
//  return 0;
//}

//int gradient_field_deform::precompute() {
//  cout << "[info] precompute for deformation...";
//  SparseMatrix<double> Ltemp = L_;
//  if ( !fixDOF_.empty() ) {
//    build_global_local_mapping((size_t)nods_.size(2), fixDOF_, g2l_);
//    rm_spmat_col_row(Ltemp, g2l_);
//  }
//  sol_.compute(Ltemp);
//  ASSERT(sol_.info() == Success);
//  return 0;
//}

int gradient_field_deform::solve_harmonic_field() { // Lf = 0
  unordered_set<size_t> nofreeDoF;
  for (auto &pid : fixDoF_)
    nofreeDoF.insert(pid);
  for (auto &pid : editDoF_) {
    nofreeDoF.insert(pid);
    hf_[pid] = 1.0;
  }
  vector<size_t> G2L;
  build_global_local_mapping((size_t)nods_.size(2), nofreeDoF, G2L);

  SparseMatrix<double> LHS = L_;
  VectorXd rhs = -LHS*hf_;
  if ( !nofreeDoF.empty() ) {
    rm_spmat_col_row(LHS, G2L);
    rm_vector_row(rhs, G2L);
  }
  SimplicialCholesky<SparseMatrix<double>> solver;
  solver.compute(LHS);
  ASSERT(solver.info() == Success);
  VectorXd df = solver.solve(rhs);
  ASSERT(solver.info() == Success);
  VectorXd Df = VectorXd::Zero(nods_.size(2));
  if ( !nofreeDoF.empty() ) {
    rc_vector_row(df, G2L, Df);
  } else {
    df = Df;
  }
  hf_ += Df;
  return 0;
}

//int gradient_field_deform::solve_for_xyz(const int xyz) {
//  matd_t x = _nods_(xyz, colon());
//  Map<VectorXd> X(&x[0], x.size());
//  VectorXd rhs;
//  calc_divergence(grad_xyz_.col(xyz), rhs);
//  rhs -= L_*X;

//  if ( !fixDOF_.empty() ) {
//    rm_vector_row(rhs, g2l_);
//  }
//  VectorXd dx = sol_.solve(rhs);
//  ASSERT(sol_.info() == Success);
//  VectorXd Dx = VectorXd::Zero(nods_.size(2));
//  if ( !fixDOF_.empty() ) {
//    rc_vector_row(dx, g2l_, Dx);
//  } else {
//    Dx = dx;
//  }
//  X += Dx;
//  _nods_(xyz, colon()) = x;

//  return 0;
//}

//int gradient_field_deform::deform(double *x) {
//  cout << "[info] deform...";
//  solve_for_xyz(0);
//  solve_for_xyz(1);
//  solve_for_xyz(2);
//  cout << "...complete!\n";
//  return 0;
//}

//int gradient_field_deform::calc_init_coord_grad() {
//  grad_xyz_ = MatrixXd::Zero(3*tris_.size(2), 3);
//  for (size_t j = 0; j < 3; ++j) {
//    matd_t coord = nods_(j, colon());
//    grad_xyz_.col(j) = G_*Map<const VectorXd>(&coord[0], coord.size());
//  }
//  return 0;
//}

//int gradient_field_deform::calc_bary_basis_grad(const mati_t &tris, const matd_t &nods, MatrixXd &gradB) {
//  gradB.setZero(3, 3*tris.size(2));
//#pragma omp parallel for
//  for (size_t i = 0; i < tris.size(2); ++i) {
//    matd_t vert = nods(colon(), tris(colon(), i));
//    calc_tri_height_vector<3>(&vert[0], &gradB(0, 3*i));
//    gradB.col(3*i+0) /= gradB.col(3*i+0).squaredNorm();
//    gradB.col(3*i+1) /= gradB.col(3*i+1).squaredNorm();
//    gradB.col(3*i+2) /= gradB.col(3*i+2).squaredNorm();
//  }
//  return 0;
//}

//int gradient_field_deform::scale_grad_fields(const double scale) {
//  grad_xyz_ *= scale;
//  return 0;
//}

//int gradient_field_deform::rotate_grad_fields(const double *axis, const double angle) {
//  Vector3d ax(axis);
//  ax /= ax.norm();
//  Matrix3d rot;
//  rot = AngleAxisd(angle, ax);
//#pragma omp parallel for
//  for (size_t i = 0; i < tris_.size(2); ++i) {
//    grad_xyz_.block<3, 3>(3*i, 0) = (rot*grad_xyz_.block<3, 3>(3*i, 0)).eval();
//  }
//  return 0;
//}

//int gradient_field_deform::reverse_grad_fields() {
//  grad_xyz_ *= -1;
//  return 0;
//}

}
