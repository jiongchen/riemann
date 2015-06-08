#include "frame_field_deform.h"

#include <iostream>
#include <fstream>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/util.h>
#include <Eigen/SVD>
#include <Eigen/UmfPackSupport>
#include <Eigen/Geometry>
#include <zjucad/matrix/itr_matrix.h>

#include "util.h"
#include "vtk.h"

using namespace std;
using namespace Eigen;
using namespace zjucad::matrix;

namespace geom_deform {

int get_ortn_basis(const double *v, double *u, double *w);

//--------------------------
//  energy definition part
//--------------------------

frame_field_deform::frame_field_deform() {}

int frame_field_deform::load_mesh(const char *file) {
  int state = 0;
  state |= jtf::mesh::load_obj(file, tris_, nods_);
  state |= build_local_frames();
  return state;
}

int frame_field_deform::save_local_frame(const char *file, const double len) const {
  mati_t lines(2, 3*tris_.size(2));
  matd_t verts(3, 4*tris_.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < tris_.size(2); ++i) {
    lines(0, 3*i+0) = 4*i+0;
    lines(1, 3*i+0) = 4*i+1;
    lines(0, 3*i+1) = 4*i+0;
    lines(1, 3*i+1) = 4*i+2;
    lines(0, 3*i+2) = 4*i+0;
    lines(1, 3*i+2) = 4*i+3;
    verts(colon(), 4*i+0) = nods_(colon(), tris_(colon(), i))*ones<double>(3, 1)/3.0;
    verts(colon(), 4*i+1) = verts(colon(), 4*i+0)+len*itr_matrix<const double*>(3, 1, &B_[i](0, 0));
    verts(colon(), 4*i+2) = verts(colon(), 4*i+0)+len*itr_matrix<const double*>(3, 1, &B_[i](0, 1));
    verts(colon(), 4*i+3) = verts(colon(), 4*i+0)+len*itr_matrix<const double*>(3, 1, &B_[i](0, 2));
  }
  ofstream os(file);
  line2vtk(os, &verts[0], verts.size(2), &lines[0], lines.size(2));
  return 0;
}

int frame_field_deform::load_constraints(const char *file) {
  ifstream is(file);
  if ( is.fail() ) {
    cerr << "[ERROR] can't load constraints\n";
    return __LINE__;
  }
  size_t fid;
  W_ = VectorXd::Zero(3*tris_.size(2));
  Vector3d v, w;
  while ( is >> fid >> v[0] >> v[1] >> v[2] >> w[0] >> w[1] >> w[2] ) {
    Matrix2d V, W;
    V(0, 0) = 0;
    V(0, 1) = 0;
    V(1, 0) = 0;
    V(1, 1) = 0;
  }
  /// init g2l
  return 0;
}

int frame_field_deform::build_local_frames() {
  matd_t normal;
  jtf::mesh::cal_face_normal(tris_, nods_, normal, true);
  B_.resize(tris_.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < tris_.size(2); ++i) {
    get_ortn_basis(&normal(0, i), &B_[i](0, 0), &B_[i](0, 1));
    std::copy(&normal(0, i), &normal(0, i)+3, &B_[i](0, 2));
  }
  return 0;
}

int frame_field_deform::save_original_mesh(const char *file) const {
  return jtf::mesh::save_obj(file, tris_, nods_);
}

int frame_field_deform::save_deformed_mesh(const char *file) const {
  return jtf::mesh::save_obj(file, tris_, _nods_);
}

int frame_field_deform::interp_frame_fields() {
  using jtf::mesh::edge2cell_adjacent;
  shared_ptr<edge2cell_adjacent> e2c(edge2cell_adjacent::create(tris_, false));

  vector<Triplet<double>> trips;
  for (auto &e : e2c->edges_) {
    pair<size_t, size_t> facet = e2c->query(e.first, e.second);
    const size_t I = facet.first;
    const size_t J = facet.second;
    // J->I
    {
      Vector3d axis = B_[J].col(2).cross(B_[I].col(2));
      axis /= axis.norm();
      double angle = surfparam::safe_acos(B_[J].col(2).dot(B_[I].col(2)));
      Matrix3d rot;
      rot = AngleAxisd(angle, axis);
      // Matrix3d R = ;
    }
    // I->J
    {
      Vector3d axis = B_[I].col(2).cross(B_[J].col(2));
      axis /= axis.norm();
    }
  }
  SparseMatrix<double> Lb;
  Lb.resize(3*tris_.size(2), 3*tris_.size(2));
  Lb.setFromTriplets(trips.begin(), trips.end());
  VectorXd b = VectorXd::Zero(3*tris_.size(2));

  surfparam::rm_spmat_col_row(Lb, g2l_);
  surfparam::rm_vector_row(b, g2l_);

  UmfPackLU<SparseMatrix<double>> sol;
  sol.compute(Lb);
  ASSERT(sol.info() == Success);
  VectorXd x = sol.solve(b);
  ASSERT(sol.info() == Success);

  surfparam::up_vector_row(x, g2l_, W_);
  return 0;
}

}
