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
#include "def.h"

using namespace std;
using namespace Eigen;
using namespace zjucad::matrix;

namespace geom_deform {

int get_ortn_basis(const double *v, double *u, double *w);

//--------------------------
//  energy definition part
//--------------------------

class deform_energy : surfparam::Functional<double>
{
public:
  deform_energy();
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, vector<Triplet<double>> *hes) const;
private:
};

frame_field_deform::frame_field_deform() {}

int frame_field_deform::load_mesh(const char *file) {
  int state = 0;
  state |= jtf::mesh::load_obj(file, tris_, nods_);
  state |= build_local_bases();
  return state;
}

int frame_field_deform::build_local_bases() {
  matd_t normal;
  jtf::mesh::cal_face_normal(tris_, nods_, normal, true);
  B_.resize(3, 3*tris_.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < tris_.size(2); ++i) {
    get_ortn_basis(&normal(0, i), &B_(0, 3*i+0), &B_(0, 3*i+1));
    std::copy(&normal(0, i), &normal(0, i)+3, &B_(0, 3*i+2));
  }
  return 0;
}

int frame_field_deform::visualize_local_bases(const char *file, const double len) const {
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
    verts(colon(), 4*i+1) = verts(colon(), 4*i+0)+len*itr_matrix<const double*>(3, 1, &B_(0, 3*i+0));
    verts(colon(), 4*i+2) = verts(colon(), 4*i+0)+len*itr_matrix<const double*>(3, 1, &B_(0, 3*i+1));
    verts(colon(), 4*i+3) = verts(colon(), 4*i+0)+len*itr_matrix<const double*>(3, 1, &B_(0, 3*i+2));
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
  W_ = VectorXd::Zero(3*tris_.size(2));     // SPD tensor field
  X_ = MatrixXd::Zero(2, 2*tris_.size(2));  // 2d cross field
  F_ = MatrixXd::Zero(3, 2*tris_.size(2));  // 3d frame field

  size_t nbr_cons;
  is >> nbr_cons;
  cout << "[INFO] constraint number: " << nbr_cons << endl;
  while ( nbr_cons-- ) {
    size_t fid;
    is >> fid;
    is >> F_(0, 2*fid+0) >> F_(1, 2*fid+0) >> F_(2, 2*fid+0);
    is >> F_(0, 2*fid+1) >> F_(1, 2*fid+1) >> F_(2, 2*fid+1);
    cons_face_.push_back(fid);
    ffc_.insert(3*fid+0);
    ffc_.insert(3*fid+1);
    ffc_.insert(3*fid+2);
  }
#pragma omp parallel for
  for (size_t i = 0; i < cons_face_.size(); ++i) {
    const size_t fid = cons_face_[i];
    Matrix2d V;
    V(0, 0) = F_.col(2*fid+0).dot(B_.col(3*fid+0));
    V(1, 0) = F_.col(2*fid+0).dot(B_.col(3*fid+1));
    V(0, 1) = F_.col(2*fid+1).dot(B_.col(3*fid+0));
    V(1, 1) = F_.col(2*fid+1).dot(B_.col(3*fid+1));
    JacobiSVD<Matrix2d> sol(V, ComputeFullU|ComputeFullV);
    Matrix2d U = sol.matrixU()*sol.matrixV().transpose();
    Matrix2d P = sol.matrixV()*sol.singularValues().asDiagonal()*sol.matrixV().transpose();
    Matrix2d W = U*P*U.transpose();
    W_[3*fid+0] = W(0, 0);
    W_[3*fid+1] = W(0, 1);
    W_[3*fid+2] = W(1, 1);
    X_.block<2, 2>(0, 2*fid) = U;
  }
  /// init global to local map
  g2l_.resize(3*tris_.size(2));
  size_t ptr = 0;
  for (size_t i = 0; i < g2l_.size(); ++i) {
    if ( ffc_.find(i) != ffc_.end() )
      g2l_[i] = -1;
    else
      g2l_[i] = ptr++;
  }
  return 0;
}

int frame_field_deform::visualize_init_frames(const char *file, const double scale) const {
  mati_t lines(2, 2*tris_.size(2));
  matd_t verts(3, 3*tris_.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < tris_.size(2); ++i) {
    lines(0, 2*i+0) = 3*i+0;
    lines(1, 2*i+0) = 3*i+1;
    lines(0, 2*i+1) = 3*i+0;
    lines(1, 2*i+1) = 3*i+2;
    verts(colon(), 3*i+0) = nods_(colon(), tris_(colon(), i))*ones<double>(3, 1)/3.0;
    verts(colon(), 3*i+1) = verts(colon(), 3*i)+scale*itr_matrix<const double*>(3, 1, &F_(0, 2*i+0));
    verts(colon(), 3*i+2) = verts(colon(), 3*i)+scale*itr_matrix<const double*>(3, 1, &F_(0, 2*i+1));
  }
  ofstream os(file);
  line2vtk(os, &verts[0], verts.size(2), &lines[0], lines.size(2));
  return 0;
}

int frame_field_deform::visualize_frame_fields(const char *file, const double scale) {
  interp_cross_fields();
  mati_t lines(2, 2*tris_.size(2));
  matd_t verts(3, 3*tris_.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < tris_.size(2); ++i) {
    lines(0, 2*i+0) = 3*i+0;
    lines(1, 2*i+0) = 3*i+1;
    lines(0, 2*i+1) = 3*i+0;
    lines(1, 2*i+1) = 3*i+2;
    Matrix2d w;
    w << W_[3*i+0], W_[3*i+1], W_[3*i+1], W_[3*i+2];
    Matrix<double, 3, 2> fr = B_.block<3, 2>(0, 3*i)* w * X_.block<2, 2>(0, 2*i);
    verts(colon(), 3*i+0) = nods_(colon(), tris_(colon(), i))*ones<double>(3, 1)/3.0;
    verts(colon(), 3*i+1) = verts(colon(), 3*i)+scale*itr_matrix<const double*>(3, 1, &fr(0, 0));
    verts(colon(), 3*i+2) = verts(colon(), 3*i)+scale*itr_matrix<const double*>(3, 1, &fr(0, 1));
  }
  ofstream os(file);
  line2vtk(os, &verts[0], verts.size(2), &lines[0], lines.size(2));
  return 0;
}

int frame_field_deform::save_original_mesh(const char *file) const {
  return jtf::mesh::save_obj(file, tris_, nods_);
}

int frame_field_deform::save_deformed_mesh(const char *file) const {
  return jtf::mesh::save_obj(file, tris_, _nods_);
}

int frame_field_deform::interp_frame_fields() {
  cout << "[INFO] interpolating tensor fields\n";
  using jtf::mesh::edge2cell_adjacent;
  shared_ptr<edge2cell_adjacent> e2c(edge2cell_adjacent::create(tris_, false));

  vector<Triplet<double>> trips;
  for (auto &e : e2c->edges_) {
    pair<size_t, size_t> facet = e2c->query(e.first, e.second);
    const size_t IJ[2] = {facet.first, facet.second};
    if ( IJ[0] == -1 || IJ[1] == -1 )
      continue;
    for (size_t k = 0; k < 1; ++k) {
      const size_t I = IJ[k];
      const size_t J = IJ[1-k];
      Vector3d axis = Vector3d(&B_(0, 3*J+2)).cross(Vector3d(&B_(0, 3*I+2)));
      axis /= axis.norm();
      double angle = surfparam::safe_acos(B_.col(3*J+2).dot(B_.col(3*I+2)));
      Matrix3d rot;
      rot = AngleAxisd(angle, axis);
      Matrix2d Rij = B_.block<3, 2>(0, 3*I).transpose()*rot*B_.block<3, 2>(0, 3*J);
      const double rija = Rij(0, 0), rijb = Rij(0, 1), rijc = Rij(1, 0), rijd = Rij(1, 1);
      trips.push_back(Triplet<double>(3*I+0, 3*I+0, -1));
      trips.push_back(Triplet<double>(3*I+0, 3*J+0, rija*rija));
      trips.push_back(Triplet<double>(3*I+0, 3*J+1, 2*rija*rijb));
      trips.push_back(Triplet<double>(3*I+0, 3*J+2, rijb*rijb));
      trips.push_back(Triplet<double>(3*I+1, 3*I+1, -1));
      trips.push_back(Triplet<double>(3*I+1, 3*J+0, rija*rijc));
      trips.push_back(Triplet<double>(3*I+1, 3*J+1, rijb*rijc+rija*rijd));
      trips.push_back(Triplet<double>(3*I+1, 3*J+2, rijb*rijd));
      trips.push_back(Triplet<double>(3*I+2, 3*I+2, -1));
      trips.push_back(Triplet<double>(3*I+2, 3*J+0, rijc*rijc));
      trips.push_back(Triplet<double>(3*I+2, 3*J+1, 2*rijc*rijd));
      trips.push_back(Triplet<double>(3*I+2, 3*J+2, rijd*rijd));
    }
  }
  const size_t dim = 3*tris_.size(2);
  SparseMatrix<double> Lb;
  Lb.resize(dim, dim);
  Lb.setFromTriplets(trips.begin(), trips.end());
  VectorXd b = -Lb * W_;

  surfparam::rm_spmat_col_row(Lb, g2l_);
  surfparam::rm_vector_row(b, g2l_);

  UmfPackLU<SparseMatrix<double>> sol;
  sol.compute(Lb);
  ASSERT(sol.info() == Success);
  VectorXd dx = sol.solve(b);
  ASSERT(sol.info() == Success);

  VectorXd DX = VectorXd::Zero(dim);
  surfparam::up_vector_row(dx, g2l_, DX);
  W_ += DX;

  return 0;
}

int frame_field_deform::visualize_tensor_fields(const char *file) {
  ofstream os(file);
  if ( os.fail() )
    return __LINE__;

  vector<string> dat_name{"spd_a", "spd_b", "spd_c", "spd_u", "spd_v"};
  MatrixXd cell_da(tris_.size(2), 5);
#pragma omp parallel for
  for (size_t i = 0; i < cell_da.rows(); ++i) {
    cell_da(i, 0) = W_[3*i+0];
    cell_da(i, 1) = W_[3*i+1];
    cell_da(i, 2) = W_[3*i+2];
    Matrix2d w;
    w << W_[3*i+0], W_[3*i+1], W_[3*i+1], W_[3*i+2];
    JacobiSVD<Matrix2d> svd(w);
    cell_da(i, 3) = svd.singularValues()[0];
    cell_da(i, 4) = svd.singularValues()[1];
  }

  os.precision(15);
  tri2vtk(os, &nods_[0], nods_.size(2), &tris_[0], tris_.size(2));

  if ( cell_da.cols() == 0 )
    return 0;
  cell_data(os, &cell_da(0, 0), cell_da.rows(), dat_name[0].c_str(),dat_name[0].c_str());
  for (size_t j = 1; j < cell_da.cols(); ++j)
    vtk_data(os, &cell_da(0, j), cell_da.rows(), dat_name[j].c_str(), dat_name[j].c_str());
  return 0;
}

bool frame_field_deform::check_spd_tensor_fields() const {
  for (size_t i = 0; i < tris_.size(2); ++i) {
    if ( W_[3*i+0] <= 0 || W_[3*i+0]*W_[3*i+2]-W_[3*i+1]*W_[3*i+1] <= 0 ) {
      return false;
    }
  }
  return true;
}

int frame_field_deform::interp_cross_fields() {
  return 0;
}

int frame_field_deform::precompute() {
  return 0;
}

int frame_field_deform::deform() {
  return 0;
}

}
