#include "frame_field_deform.h"

#include <iostream>
#include <fstream>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/util.h>
#include <Eigen/SVD>
#include <Eigen/UmfPackSupport>
#include <Eigen/Geometry>
#include <zjucad/matrix/itr_matrix.h>

#include "def.h"
#include "util.h"
#include "vtk.h"
#include "cotmatrix.h"

using namespace std;
using namespace Eigen;
using namespace zjucad::matrix;
using namespace jtf::mesh;

namespace riemann {

//--------------------------
//- energy definition part -
//--------------------------

class deform_energy : public riemann::Functional<double>
{
public:
  typedef zjucad::matrix::matrix<size_t> mati_t;
  typedef zjucad::matrix::matrix<double> matd_t;
  deform_energy(const mati_t &tris, const matd_t &nods, const vector<Matrix3d> &Winv, const double w)
    : tris_(tris), nods_(nods), Winv_(Winv), w_(w) {
    // init optimal rotation
    Q_.resize(tris_.size(2));
#pragma omp parallel for
    for (size_t i = 0; i < Q_.size(); ++i)
      Q_[i] = Matrix3d::Identity();

    face_cot_.resize(3, tris_.size(2));
#pragma omp parallel for
    for (size_t i = 0; i < face_cot_.size(2); ++i) {
      matd_t vert = nods_(colon(), tris_(colon(), i));
      face_cot_(0, i) = riemann::cal_cot_val(&vert(0, 0), &vert(0, 2), &vert(0, 1));
      face_cot_(1, i) = riemann::cal_cot_val(&vert(0, 1), &vert(0, 0), &vert(0, 2));
      face_cot_(2, i) = riemann::cal_cot_val(&vert(0, 2), &vert(0, 1), &vert(0, 0));
    }

    e2c_.reset(edge2cell_adjacent::create(tris_, false));
    cotval_.resize(2, e2c_->edges_.size());
#pragma omp parallel for
    for (size_t i = 0; i < e2c_->edges_.size(); ++i) {
      const size_t ei = e2c_->edges_[i].first;
      const size_t ej = e2c_->edges_[i].second;
      pair<size_t, size_t> fa = e2c_->query(ei, ej);
      const size_t face[2] = {fa.first, fa.second};
      for (size_t k = 0; k < 2; ++k) {
        if ( face[k] == -1 ) {
          cotval_(k, i) = 0.0;
          continue;
        }
        const size_t eo = sum(tris_(colon(), face[k]))-ei-ej;
        cotval_(k, i) = riemann::cal_cot_val(&nods_(0, ei), &nods_(0, eo), &nods_(0, ej));
      }
    }
  }
  size_t Nx() const {
    return nods_.size();
  }
  int Val(const double *x, double *val) const {
    Map<const MatrixXd> X(x, 3, Nx()/3);
    Map<const MatrixXd> P(&nods_[0], 3, Nx()/3);
    for (size_t i = 0; i < e2c_->edges_.size(); ++i) {
      const size_t ei = e2c_->edges_[i].first;
      const size_t ej = e2c_->edges_[i].second;
      pair<size_t, size_t> fa = e2c_->query(ei, ej);
      const size_t face[2] = {fa.first, fa.second};
      for (size_t k = 0; k < 2; ++k) {
        if ( face[k] == -1 )
          continue;
        *val += w_*cotval_(k, i)*
            (X.col(ei)-X.col(ej)-Q_[face[k]]*Winv_[face[k]]*(P.col(ei)-P.col(ej))).squaredNorm();
      }
    }
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    Map<const MatrixXd> X(x, 3, Nx()/3);
    Map<const MatrixXd> P(&nods_[0], 3, Nx()/3);
    Map<MatrixXd> grad(gra, 3, Nx()/3);
    for (size_t i = 0; i < e2c_->edges_.size(); ++i) {
      const size_t ei = e2c_->edges_[i].first;
      const size_t ej = e2c_->edges_[i].second;
      pair<size_t, size_t> fa = e2c_->query(ei, ej);
      const size_t face[2] = {fa.first, fa.second};
      for (size_t k = 0; k < 2; ++k) {
        if ( face[k] == -1 )
          continue;
        Vector3d temp = 2.0*w_*cotval_(k, i)*
            (X.col(ei)-X.col(ej)-Q_[face[k]]*Winv_[face[k]]*(P.col(ei)-P.col(ej)));
        grad.col(ei) += temp;
        grad.col(ej) -= temp;
      }
    }
    return 0;
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    for (size_t i = 0; i < e2c_->edges_.size(); ++i) {
      const size_t ei = e2c_->edges_[i].first;
      const size_t ej = e2c_->edges_[i].second;
      pair<size_t, size_t> fa = e2c_->query(ei, ej);
      const size_t face[2] = {fa.first, fa.second};
      for (size_t k = 0; k < 2; ++k) {
        if ( face[k] == -1 )
          continue;
        const double wgt = 2.0*w_*cotval_(k, i);
        riemann::add_diag_block(ei, ei, +wgt, hes);
        riemann::add_diag_block(ei, ej, -wgt, hes);
        riemann::add_diag_block(ej, ej, +wgt, hes);
        riemann::add_diag_block(ej, ei, -wgt, hes);
      }
    }
    return 0;
  }
  int EvaluateOptimalRotation(const double *x) {
    Map<const MatrixXd> X(x, 3, Nx()/3);
    Map<const MatrixXd> P(&nods_[0], 3, Nx()/3);
#pragma omp parallel for
    for (size_t i = 0; i < tris_.size(2); ++i) {
      Matrix3d curr_pts, rest_pts;
      Vector3d wgt(face_cot_(0, i), face_cot_(1, i), face_cot_(2, i));

      rest_pts.col(0) = P.col(tris_(1, i))-P.col(tris_(0, i));
      rest_pts.col(1) = P.col(tris_(2, i))-P.col(tris_(1, i));
      rest_pts.col(2) = P.col(tris_(0, i))-P.col(tris_(2, i));
      rest_pts = (Winv_[i]*rest_pts).eval();

      curr_pts.col(0) = X.col(tris_(1, i))-X.col(tris_(0, i));
      curr_pts.col(1) = X.col(tris_(2, i))-X.col(tris_(1, i));
      curr_pts.col(2) = X.col(tris_(0, i))-X.col(tris_(2, i));

      Matrix3d S = curr_pts*wgt.asDiagonal()*rest_pts.transpose();
      JacobiSVD<Matrix3d> svd(S, ComputeFullU|ComputeFullV);
      Matrix3d su = svd.matrixU();
      Matrix3d sv = svd.matrixV();
      Q_[i] = su*sv.transpose();
      if ( Q_[i].determinant() < 0 ) {
        su.col(2) = -su.col(2);
        Q_[i] = su*sv.transpose();
      }
    }
    return 0;
  }
  int EvaluateOptimalRotationOther(const double *x) {
    Map<const MatrixXd> X(x, 3, Nx()/3);
    Map<const MatrixXd> P(&nods_[0], 3, Nx()/3);
#pragma omp parallel for
    for (size_t i = 0; i < tris_.size(2); ++i) {
      Matrix3d curr_pts, rest_pts;
      Vector3d wgt(face_cot_(0, i), face_cot_(1, i), face_cot_(2, i));
      Vector3d curr_cent, rest_cent;

      rest_pts.col(0) = P.col(tris_(1, i))-P.col(tris_(0, i));
      rest_pts.col(1) = P.col(tris_(2, i))-P.col(tris_(1, i));
      rest_pts.col(2) = P.col(tris_(0, i))-P.col(tris_(2, i));
      rest_pts = (Winv_[i]*rest_pts).eval();
      rest_cent = rest_pts*wgt/wgt.sum();

      curr_pts.col(0) = X.col(tris_(1, i))-X.col(tris_(0, i));
      curr_pts.col(1) = X.col(tris_(2, i))-X.col(tris_(1, i));
      curr_pts.col(2) = X.col(tris_(0, i))-X.col(tris_(2, i));
      curr_cent = curr_pts*wgt/wgt.sum();

      Matrix3d XX = rest_pts-rest_cent*Vector3d::Ones().transpose();
      Matrix3d YY = curr_pts-curr_cent*Vector3d::Ones().transpose();
      Matrix3d SS = XX*wgt.asDiagonal()*YY.transpose();
      JacobiSVD<Matrix3d> svd(SS, ComputeFullU|ComputeFullV);
      Matrix3d su = svd.matrixU();
      Matrix3d sv = svd.matrixV();
      Matrix3d I = Matrix3d::Identity();
      I(2, 2) = (sv*su.transpose()).determinant();
      Q_[i] = sv*I*su.transpose();
    }
    return 0;
  }
private:
  const mati_t &tris_;
  const matd_t &nods_;
  const vector<Matrix3d> Winv_;
  const double w_;
  vector<Matrix3d> Q_;
  matd_t cotval_;
  matd_t face_cot_;
  shared_ptr<edge2cell_adjacent> e2c_;
};

class smooth_energy : public riemann::Functional<double>
{
public:
  typedef zjucad::matrix::matrix<size_t> mati_t;
  typedef zjucad::matrix::matrix<double> matd_t;
  smooth_energy(const mati_t &tris, const matd_t &nods, const double w)
    : tris_(tris), nods_(nods), w_(w) {
    riemann::cotmatrix(tris_, nods_, 3, &L_);
    LtL_ = L_.transpose()*L_;
  }
  size_t Nx() const {
    return nods_.size();
  }
  int Val(const double *x, double *val) const {
    Map<const VectorXd> X(x, Nx());
    Map<const VectorXd> p(&nods_[0], Nx());
    *val += w_*(L_*(X-p)).squaredNorm();
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    Map<const VectorXd> X(x, Nx());
    Map<const VectorXd> p(&nods_[0], Nx());
    Map<VectorXd> grad(gra, Nx());
    grad += 2.0*w_*LtL_*(X-p);
    return 0;
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    for (size_t j = 0; j < LtL_.outerSize(); ++j) {
      for (SparseMatrix<double>::InnerIterator it(LtL_, j); it; ++it)
        hes->push_back(Triplet<double>(it.row(), it.col(), 2*w_*it.value()));
    }
    return 0;
  }
private:
  const mati_t &tris_;
  const matd_t &nods_;
  const double w_;
  SparseMatrix<double> L_, LtL_;
};

frame_field_deform::frame_field_deform()
  : max_iter_(20000),
    tolerance_(1e-12),
    lambda_(0.1),
    perturb_(0.1) {}

frame_field_deform::frame_field_deform(const boost::property_tree::ptree &pt) {
  max_iter_  =  pt.get<size_t>("max_iter");
  tolerance_ =  pt.get<double>("tolerance");
  lambda_    =  pt.get<double>("lambda");
  perturb_   =  pt.get<double>("perturb");
}

int frame_field_deform::load_mesh(const char *file) {
  int state = 0;
  state |= jtf::mesh::load_obj(file, tris_, nods_);
  state |= build_local_bases(tris_, nods_, B_);
  return state;
}

int frame_field_deform::build_local_bases(const mati_t &tris, const matd_t &nods, MatrixXd &B) {
  matd_t normal;
  jtf::mesh::cal_face_normal(tris, nods, normal, true);
  B.resize(3, 3*tris.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < tris.size(2); ++i) {
    matd_t u = nods(colon(), tris(1, i))-nods(colon(), tris(0, i));
    u /= norm(u);
    matd_t v = cross(normal(colon(), i), u);
    v /= norm(v);
    std::copy(&u[0], &u[0]+3, &B(0, 3*i+0));
    std::copy(&v[0], &v[0]+3, &B(0, 3*i+1));
    std::copy(&normal(0, i), &normal(0, i)+3, &B(0, 3*i+2));
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
  os.close();
  return 0;
}

int frame_field_deform::load_constraints(const char *file) {
  ifstream is(file);
  if ( is.fail() ) {
    cerr << "[error] can't load constraints\n";
    return __LINE__;
  }
  W_ = VectorXd::Zero(3*tris_.size(2));     // SPD tensor field
  X_ = MatrixXd::Zero(2, 2*tris_.size(2));  // 2d cross field
  F_ = MatrixXd::Zero(3, 2*tris_.size(2));  // 3d frame field

  size_t nbr_cons;
  is >> nbr_cons;
  cout << "[info] constraint number: " << nbr_cons << endl;
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
  is.close();
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
  os.close();
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
  os.close();
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
    for (size_t k = 0; k < 2; ++k) {
      const size_t I = IJ[k];
      const size_t J = IJ[1-k];
      Matrix3d rot;
      double angle = riemann::safe_acos(B_.col(3*J+2).dot(B_.col(3*I+2)));
      if ( std::fabs(angle) < 1e-16 ) {
        rot = Matrix3d::Identity();
      } else if ( std::fabs(angle - 3.141592653589793238462643383279502884) < 1e-16 ) {
        rot = -Matrix3d::Identity();
      } else {
        Vector3d axis = Vector3d(&B_(0, 3*J+2)).cross(Vector3d(&B_(0, 3*I+2)));
        axis /= axis.norm();
        rot = AngleAxisd(angle, axis);
      }
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

  riemann::rm_spmat_col_row(Lb, g2l_);
  riemann::rm_vector_row(b, g2l_);

  UmfPackLU<SparseMatrix<double>> sol;
  sol.compute(Lb);
  ASSERT(sol.info() == Success);
  VectorXd dx = sol.solve(b);
  ASSERT(sol.info() == Success);

  VectorXd DX = VectorXd::Zero(dim);
  riemann::rc_vector_row(dx, g2l_, DX);
  W_ += DX;

  return 0;
}

int frame_field_deform::visualize_tensor_fields(const char *file) {
  ofstream os(file);
  if ( os.fail() ) {
    cerr << "[error] can't open write " << file << endl;
    return __LINE__;
  }

  vector<string> dat_name{"spd_a", "spd_b", "spd_c", "spd_u", "spd_v"};
  MatrixXd cell_dat(tris_.size(2), 5);
#pragma omp parallel for
  for (size_t i = 0; i < cell_dat.rows(); ++i) {
    cell_dat(i, 0) = W_[3*i+0];
    cell_dat(i, 1) = W_[3*i+1];
    cell_dat(i, 2) = W_[3*i+2];
    Matrix2d w;
    w << W_[3*i+0], W_[3*i+1], W_[3*i+1], W_[3*i+2];
    JacobiSVD<Matrix2d> svd(w);
    cell_dat(i, 3) = svd.singularValues()[0];
    cell_dat(i, 4) = svd.singularValues()[1];
  }

  os.precision(15);
  tri2vtk(os, &nods_[0], nods_.size(2), &tris_[0], tris_.size(2));
  if ( cell_dat.cols() == 0 )
    return 0;
  cell_data(os, &cell_dat(0, 0), cell_dat.rows(), dat_name[0].c_str(),dat_name[0].c_str());
  for (size_t j = 1; j < cell_dat.cols(); ++j)
    vtk_data(os, &cell_dat(0, j), cell_dat.rows(), dat_name[j].c_str(), dat_name[j].c_str());
  os.close();
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
  // build transformation w.r.t global bases
  std::vector<Matrix3d> Winv(tris_.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < Winv.size(); ++i) {
    Matrix2d w;
    w << W_[3*i+0], W_[3*i+1], W_[3*i+1], W_[3*i+2];
    FullPivLU<Matrix2d> sol(w);
    if ( sol.isInvertible() ) {
      Winv[i] = B_.block<3, 2>(0, 3*i)*sol.inverse()*B_.block<3, 2>(0, 3*i).transpose();
    } else {
      cerr << "[error] W" << i << " is not invertible\n";
      exit(EXIT_FAILURE);
    }
  }

  // construct energy
  buff_.resize(2);
  buff_[DEFORM] = std::make_shared<deform_energy>(tris_, nods_, Winv, 1.0-lambda_);
  buff_[SMOOTH] = std::make_shared<smooth_energy>(tris_, nods_, lambda_);
  try {
    e_ = std::make_shared<riemann::energy_t<double>>(buff_);
  } catch ( exception &e ) {
    cerr << "[error] " << e.what() << endl;
    exit(EXIT_FAILURE);
  }

  // query constant system matrix
  vector<Triplet<double>> trips;
  e_->Hes(nullptr, &trips);
  LHS_.resize(e_->Nx(), e_->Nx());
  LHS_.setFromTriplets(trips.begin(), trips.end());
  sol_.compute(LHS_);
  ASSERT(sol_.info() == Success);

  return 0;
}

int frame_field_deform::deform() {
  _nods_ = nods_;
#pragma omp parallel for
  for (size_t j = 0; j < _nods_.size(2); ++j)
    for (size_t i = 0; i < _nods_.size(1); ++i)
      _nods_(i, j) += (double(rand())/double(RAND_MAX))*10e-4*perturb_;

  Map<VectorXd> X(&_nods_[0], _nods_.size());
  for (size_t iter = 0; iter < max_iter_; ++iter) {
    // query energy value
    if ( iter % 100 == 0 ) {
      cout << "[info] iteration " << iter << endl;
      double val = 0;
      e_->Val(&X[0], &val);
      cout << "[info] energy: " << val << "\n";
    }

    // assemble rhs
    VectorXd rhs = VectorXd::Zero(e_->Nx());
    e_->Gra(&X[0], &rhs[0]);
    rhs = -rhs;
    if ( iter % 100 == 0 )
      cout << "[info] gradient norm: " << rhs.norm() << "\n\n";

    VectorXd dx = sol_.solve(rhs);
    ASSERT(sol_.info() == Success);

    double x0norm = X.norm();
    X += dx;
    std::dynamic_pointer_cast<deform_energy>(buff_[DEFORM])
        ->EvaluateOptimalRotation(&X[0]);

    // convergence test
    if ( dx.norm() <= tolerance_*x0norm ) {
      cout << "[info] converged!\n";
      break;
    }
  }
  return 0;
}

int frame_field_deform::calc_defo_grad_oper() {
  D_.resize(3, 3*tris_.size(2));
  matd_t normal;
  jtf::mesh::cal_face_normal(tris_, nods_, normal, true);
#pragma omp parallel for
  for (size_t i = 0; i < tris_.size(2); ++i) {
    matd_t D(3, 3);
    D(colon(), colon(0, 1)) = nods_(colon(), tris_(colon(1, 2), i))-nods_(colon(), tris_(colon(0, 1), i));
    D(colon(), 2) = normal(colon(), i);
    std::copy(D.begin(), D.end(), &D_(0, 3*i));
    FullPivLU<Matrix3d> flu(D_.block<3, 3>(0, 3*i));
    if ( flu.isInvertible() )
      D_.block<3, 3>(0, 3*i) = flu.inverse();
    else
      cerr << "[info] [e0, e1, n] is not invertible\n";
  }
  return 0;
}

int frame_field_deform::gen_cross_field() {
  X_ = MatrixXd::Zero(3, 2*tris_.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < cons_face_.size(); ++i) {
    const size_t fid = cons_face_[i];
    matd_t temp = zeros<double>(3, 3);
    temp(colon(), colon(0, 1)) = _nods_(colon(), tris_(colon(1, 2), fid))-_nods_(colon(), tris_(colon(0, 1), fid));
    Matrix3d DG = Matrix3d(temp.begin())*D_.block<3, 3>(0, 3*fid);
    X_.col(2*fid+0) = DG*F_.col(2*fid+0);
    X_.col(2*fid+1) = DG*F_.col(2*fid+1);
  }
  return 0;
}

int frame_field_deform::visualize_cross_fields(const char *file, const double scale) const {
  mati_t lines(2, 2*tris_.size(2));
  matd_t verts(3, 3*tris_.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < tris_.size(2); ++i) {
    lines(0, 2*i+0) = 3*i+0;
    lines(1, 2*i+0) = 3*i+1;
    lines(0, 2*i+1) = 3*i+0;
    lines(1, 2*i+1) = 3*i+2;
    verts(colon(), 3*i+0) = _nods_(colon(), tris_(colon(), i))*ones<double>(3, 1)/3.0;
    verts(colon(), 3*i+1) = verts(colon(), 3*i)+scale*itr_matrix<const double*>(3, 1, &X_(0, 2*i+0));
    verts(colon(), 3*i+2) = verts(colon(), 3*i)+scale*itr_matrix<const double*>(3, 1, &X_(0, 2*i+1));
  }
  ofstream os(file);
  line2vtk(os, &verts[0], verts.size(2), &lines[0], lines.size(2));
  os.close();
  return 0;
}

}
