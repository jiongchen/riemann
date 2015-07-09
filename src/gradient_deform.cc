#include "gradient_deform.h"

#include <zjucad/matrix/itr_matrix.h>
#include <jtflib/mesh/io.h>

#include "cotmatrix.h"
#include "grad_operator.h"
#include "vtk.h"

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
  // compute gradient operator
  calc_grad_operator(tris_, nods_, &G_);
  // compute initial coordinate gradient
  calc_init_coord_grad();
  // compute the gradient of barycentric basis function
  calc_bary_basis_grad();
  // comupte triangle areas
  calc_element_area();
  return 0;
}

int gradient_field_deform::see_coord_grad_fields(const char *filename, const int xyz) const {
  mati_t line(2, tris_.size(2));
  matd_t vert(3, 2*tris_.size(2));
  line(0, colon()) = colon(0, tris_.size(2)-1);
  line(1, colon()) = colon(tris_.size(2), 2*tris_.size(2)-1);
  itr_matrix<const double *> g(3, tris_.size(2), &grad_xyz_(0, xyz));
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
    for (size_t j = 0; j < 3; ++i) {
      div[tris_(j, i)] += gradB_.col(3*i+j).dot(vf.segment<3>(3*i))*area_[i];
    }
  }
  return 0;
}

}
