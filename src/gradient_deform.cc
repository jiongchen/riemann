#include "gradient_deform.h"

#include <zjucad/matrix/itr_matrix.h>
#include <jtflib/mesh/io.h>

#include "cotmatrix.h"
#include "grad_operator.h"
#include "vtk.h"
#include "util.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

namespace geom_deform {

gradient_field_deform::gradient_field_deform() {

}

int gradient_field_deform::load_origin_model(const char *filename) {
  int rtn = jtf::mesh::load_obj(filename, tris_, nods_);
  return rtn;
}

int gradient_field_deform::save_origin_model(const char *filename) const {
  int rtn = jtf::mesh::save_obj(filename, tris_, nods_);
  return rtn;
}

int gradient_field_deform::save_deformed_model(const char *filename) const {
  int rtn = jtf::mesh::save_obj(filename, tris_, _nods_);
  return rtn;
}

int gradient_field_deform::calc_init_coord_grad() {
  grad_xyz_ = MatrixXd::Zero(3*tris_.size(2), 3);
  for (size_t j = 0; j < 3; ++j) {
    matd_t coord = nods_(j, colon());
    grad_xyz_.col(j) = G_*Map<const VectorXd>(&coord[0], coord.size());
  }
  return 0;
}

int gradient_field_deform::calc_bary_basis_grad() {
  gradB_ = MatrixXd::Zero(3, 3*tris_.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < tris_.size(2); ++i) {
    matd_t vert = nods_(colon(), tris_(colon(), i));
    calc_tri_height_vector(&vert[0], &gradB_(0, 3*i));
    gradB_.col(3*i+0) /= gradB_.col(3*i+0).squaredNorm();
    gradB_.col(3*i+1) /= gradB_.col(3*i+1).squaredNorm();
    gradB_.col(3*i+2) /= gradB_.col(3*i+2).squaredNorm();
  }
  return 0;
}

int gradient_field_deform::calc_element_area() {
  area_.resize(tris_.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < tris_.size(2); ++i) {
    matd_t edge = nods_(colon(), tris_(colon(1, 2), i))-nods_(colon(), tris_(colon(0, 1), i));
    area_[i] = norm(cross(edge(colon(), 0), edge(colon(), 1)))/2.0;
  }
  return 0;
}

int gradient_field_deform::init() {
  // compute Laplacian operator
  surfparam::cotmatrix(tris_, nods_, 1, &L_);
  L_ *= -1.0;
  // compute gradient operator
  calc_grad_operator(tris_, nods_, &G_);
  // compute initial coordinate gradient
  calc_init_coord_grad();
  // compute the gradient of barycentric basis function
  calc_bary_basis_grad();
  // comupte triangle areas
  calc_element_area();
  // assign deform vertices
  _nods_ = nods_;
  // initialize harmonic scalar field
  hf_ = VectorXd::Zero(nods_.size(2));
  return 0;
}

int gradient_field_deform::see_coord_grad_fields(const char *filename, const int xyz) const {
  mati_t line(2, tris_.size(2));
  matd_t vert(3, 2*tris_.size(2));
  line(0, colon()) = colon(0, tris_.size(2)-1);
  line(1, colon()) = colon(tris_.size(2), 2*tris_.size(2)-1);
  itr_matrix<const double *> g(3, grad_xyz_.rows()/3, &grad_xyz_(0, xyz));
#pragma omp parallel for
  for (size_t i = 0; i < tris_.size(2); ++i) {
    vert(colon(), i) = nods_(colon(), tris_(colon(), i))*ones<double>(3, 1)/3.0;
    vert(colon(), i+tris_.size(2)) = vert(colon(), i)+g(colon(), i);
  }
  ofstream os(filename);
  line2vtk(os, &vert[0], vert.size(2), &line[0], line.size(2));
  return 0;
}

int gradient_field_deform::calc_divergence(const VectorXd &vf, VectorXd &div) const {
  div = VectorXd::Zero(nods_.size(2));
  for (size_t i = 0; i < tris_.size(2); ++i) {
    for (size_t j = 0; j < 3; ++j) {
      div[tris_(j, i)] += gradB_.col(3*i+j).dot(vf.segment<3>(3*i))*area_[i];
    }
  }
  return 0;
}

int gradient_field_deform::scale_grad_fields(const double scale) {
  grad_xyz_ *= scale;
  return 0;
}

int gradient_field_deform::set_fixed_verts(const vector<size_t> &idx) {
  for (auto &id : idx) {
    fix_dofs_.insert(id);
  }
  return 0;
}

int gradient_field_deform::manipualte(const size_t idx, const double *u) {
  fix_dofs_.insert(idx);
  hf_[idx] = 1.0;
  _nods_(colon(), idx) += itr_matrix<const double*>(3, 1, u);
  return 0;
}

int gradient_field_deform::precompute() {
  cout << "[info] precompute for deformation...";
  SparseMatrix<double> Ltemp = L_;
  if ( !fix_dofs_.empty() ) {
    surfparam::build_global_local_mapping((size_t)nods_.size(2), fix_dofs_, g2l_);
    surfparam::rm_spmat_col_row(Ltemp, g2l_);
  }
  sol_.compute(Ltemp);
  ASSERT(sol_.info() == Success);

  // solve for harmonic fields: L(f0+df)=0
  VectorXd RHS = -L_*hf_;
  if ( !fix_dofs_.empty() )
    surfparam::rm_vector_row(RHS, g2l_);
  VectorXd df = sol_.solve(RHS);
  ASSERT(sol_.info() == Success);
  VectorXd Df = VectorXd::Zero(nods_.size(2));
  if ( !fix_dofs_.empty() )
    surfparam::rc_vector_row(df, g2l_, Df);
  else
    df = Df;
  hf_ += Df;

  // apply blended local deformation
  // .....

  cout << "...complete!\n";
  return 0;
}

int gradient_field_deform::see_harmonic_field(const char *filename) const {
  ofstream os(filename);
  if ( os.fail() )
    return __LINE__;
  tri2vtk(os, &nods_[0], nods_.size(2), &tris_[0], tris_.size(2));
  point_data(os, hf_.data(), hf_.rows(), "hf", "hf");
  os.close();
  return 0;
}

int gradient_field_deform::solve_for_xyz(const int xyz) {
  matd_t x = _nods_(xyz, colon());
  Map<VectorXd> X(&x[0], x.size());
  VectorXd rhs;
  calc_divergence(grad_xyz_.col(xyz), rhs);
  rhs -= L_*X;

  if ( !fix_dofs_.empty() ) {
    surfparam::rm_vector_row(rhs, g2l_);
  }
  VectorXd dx = sol_.solve(rhs);
  ASSERT(sol_.info() == Success);
  VectorXd Dx = VectorXd::Zero(nods_.size(2));
  if ( !fix_dofs_.empty() ) {
    surfparam::rc_vector_row(dx, g2l_, Dx);
  } else {
    Dx = dx;
  }
  X += Dx;
  _nods_(xyz, colon()) = x;

  return 0;
}

int gradient_field_deform::deform() {
  cout << "[info] deform...";
  solve_for_xyz(0);
  solve_for_xyz(1);
  solve_for_xyz(2);
  cout << "...complete!\n";
  return 0;
}

}
