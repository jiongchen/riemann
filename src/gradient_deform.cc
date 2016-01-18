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

static void calc_shape_func_grad(const mati_t &tris, const matd_t &nods, MatrixXd &gradS) {
  gradS.setZero(9, tris.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < tris.size(2); ++i) {
    matd_t vert = nods(colon(), tris(colon(), i));          // [v0, v1, v2]
    calc_tri_linear_basis_grad<3>(&vert[0], &gradS(0, i));  // [g1, g2, g3]
  }
}

static void calc_surf_scalar_field_grad(const mati_t &tris, const matd_t &u, const MatrixXd &gradS, double *gradU) {
  ASSERT(u.size(1) == 1 || u.size(2) == 1);
  Map<VectorXd> gU(gradU, 3*tris.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < tris.size(2); ++i) {
    matd_t ui = u(tris(colon(), i));
    gU.segment<3>(3*i) = Map<const Matrix3d>(&gradS(0, i))*Map<const Vector3d>(&ui[0]);
  }
}

gradient_field_deform::gradient_field_deform(const mati_t &tris, const matd_t &nods)
  : dim_(nods.size()), tris_(tris), nods_(nods), ID_(Quaterniond::Identity()) {
  // laplacian
  cotmatrix(tris_, nods_, 1, &L_);
  L_ *= -1.0;
  // gradient of shape function
  calc_shape_func_grad(tris_, nods_, gradShape_);
  // initial gradient field
  Gxyz_.setZero(3*tris_.size(2), nods_.size(2)); {
    matd_t x = nods_(0, colon());
    calc_surf_scalar_field_grad(tris_, x, gradShape_, &Gxyz_(0, 0));
    matd_t y = nods_(1, colon());
    calc_surf_scalar_field_grad(tris_, y, gradShape_, &Gxyz_(0, 1));
    matd_t z = nods_(2, colon());
    calc_surf_scalar_field_grad(tris_, z, gradShape_, &Gxyz_(0, 2));
  }
  // area
  calc_tris_area(tris_, nods_, area_);
  // harmonic field
  hf_.setZero(nods_.size(2));
  // transform field
  vs_.resize(nods_.size(2));
  vR_.resize(nods_.size(2));
  fs_.resize(tris_.size(2));
  fR_.resize(tris_.size(2));
}

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

void gradient_field_deform::prescribe_uniform_transform(const Quaterniond &q, const double s) {
  q_ = q;
  s_ = s;
}

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

void gradient_field_deform::propogate_transform() {
#pragma omp parallel for
  for (size_t i = 0; i < nods_.size(2); ++i) {
    vR_[i] = ID_.slerp(hf_[i], q_);
    vs_[i] = hf_[i]*s_+(1.0-hf_[i])*1.0;
  }
}

void gradient_field_deform::interp_face_transform() {
#pragma omp parallel for
  for (size_t i = 0; i < tris_.size(2); ++i) {
    Quaterniond q12 = vR_[tris_(0, i)].slerp(0.5, vR_[tris_(1, i)]);
    fR_[i] = q12.slerp(1.0/3, vR_[tris_(2, i)]);
    fs_[i] = (vs_[tris_(0, i)]+vs_[tris_(1, i)]+vs_[tris_(2, i)])/3.0;
  }
}

void gradient_field_deform::rotate_surf_piecewise(mati_t &rtris, matd_t &rnods) const {
  rtris.resize(3, tris_.size(2));
  rnods.resize(3, 3*tris_.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < rtris.size(2); ++i) {
    Matrix3d rot = fR_[i].toRotationMatrix();
    itr_matrix<const double *> R(3, 3, rot.data());
    matd_t vert = nods_(colon(), tris_(colon(), i));
    matd_t bc = vert*ones<double>(3, 1)/3.0;
    rtris(colon(), i) = colon(3*i, 3*i+2);
    rnods(colon(), 3*i+0) = bc+R*(nods_(colon(), tris_(0, i))-bc);
    rnods(colon(), 3*i+1) = bc+R*(nods_(colon(), tris_(1, i))-bc);
    rnods(colon(), 3*i+2) = bc+R*(nods_(colon(), tris_(2, i))-bc);
  }
}

void gradient_field_deform::update_gradient_field() {
#pragma omp parallel for
  for (size_t i = 0; i < tris_.size(2); ++i) {
    Matrix3d rot = fR_[i].toRotationMatrix();
    Gxyz_.block<3, 3>(3*i, 0) = (rot.transpose()*Gxyz_.block<3, 3>(3*i, 0)).eval();
  }
}

void gradient_field_deform::prefactorize() {
  size_t cnt = 0;
  g2l_.resize(nods_.size(2));
  for (size_t i = 0; i < g2l_.size(); ++i) {
    if ( fixDoF_.find(i) != fixDoF_.end() )
      g2l_[i] = -1;
    else
      g2l_[i] = cnt++;
  }
  SparseMatrix<double> LHS = L_;
  rm_spmat_col_row(LHS, g2l_);
  solver_.setMode(SimplicialCholeskyLLT);
  solver_.compute(LHS);
  ASSERT(solver_.info() == Success);
}

int gradient_field_deform::precompute() {
  solve_harmonic_field();     // update hf
  propogate_transform();      // update vs and vR
  interp_face_transform();    // update fs and fR
//  rotate_surf_piecewise();
  update_gradient_field();    // update Gxyz
  prefactorize();             // update solver
  return 0;
}

int gradient_field_deform::calc_div(const VectorXd &vf, VectorXd &div) const {
  div.setZero(nods_.size(2));
  for (size_t i = 0; i < tris_.size(2); ++i) {
    div[tris_(0, i)] += gradShape_.block<3, 1>(0, i).dot(vf.segment<3>(3*i))*area_[i];
    div[tris_(1, i)] += gradShape_.block<3, 1>(3, i).dot(vf.segment<3>(3*i))*area_[i];
    div[tris_(2, i)] += gradShape_.block<3, 1>(6, i).dot(vf.segment<3>(3*i))*area_[i];
  }
  return 0;
}

int gradient_field_deform::solve_xyz(double *x, const int xyz) const {
  Map<MatrixXd> X(x, 3, nods_.size(2));
  VectorXd xstar = X.row(xyz).transpose(), rhs;

  calc_div(Gxyz_.col(xyz), rhs);
  rhs -= L_*xstar;

  if ( !fixDoF_.empty() ) {
    rm_vector_row(rhs, g2l_);
  }
  VectorXd dx = solver_.solve(rhs);
  ASSERT(solver_.info() == Success);
  VectorXd Dx = VectorXd::Zero(nods_.size(2));
  if ( !fixDoF_.empty() ) {
    rc_vector_row(dx, g2l_, Dx);
  } else {
    Dx = dx;
  }
  xstar += Dx;

  X.row(xyz) = xstar.transpose();
  return 0;
}

int gradient_field_deform::deform(double *x) const {
  cout << "[info] deform...";
  solve_xyz(x, 0);
  solve_xyz(x, 1);
  solve_xyz(x, 2);
  cout << "...complete!\n";
  return 0;
}

}
