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
#include "cotmatrix.h"

using namespace std;
using namespace Eigen;
using namespace zjucad::matrix;
using namespace jtf::mesh;

namespace geom_deform {

//--------------------------
//- energy definition part -
//--------------------------

class deform_energy : public surfparam::Functional<double>
{
public:
  typedef zjucad::matrix::matrix<size_t> mati_t;
  typedef zjucad::matrix::matrix<double> matd_t;
  deform_energy(const mati_t &tris, const matd_t &nods, const vector<Matrix3d> &Winv, const double w)
    : tris_(tris), nods_(nods), Winv_(Winv), w_(w) {
    // init Q
    Q_.resize(tris_.size(2));
#pragma omp parallel for
    for (size_t i = 0; i < Q_.size(); ++i)
      Q_[i] = Matrix3d::Identity();

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
        cotval_(k, i) = surfparam::cal_cot_val(&nods_(0, ei), &nods_(0, eo), &nods_(0, ej));
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
        surfparam::add_diag_block(ei, ei, +wgt, hes);
        surfparam::add_diag_block(ei, ej, -wgt, hes);
        surfparam::add_diag_block(ej, ej, +wgt, hes);
        surfparam::add_diag_block(ej, ei, -wgt, hes);
      }
    }
    return 0;
  }
  int UpdateRotation(const double *x) {
    itr_matrix<const double *> X(3, Nx()/3, x);
    static const matrix<double> op = identity_matrix<double>(3)- ones<double>(3, 3) / 3.0;
#pragma omp parallel for
    for (size_t i = 0; i < tris_.size(2); ++i) {
      matd_t xx = X(colon(), tris_(colon(), i))*op;
      matd_t yy = nods_(colon(), tris_(colon(), i))*op;
      matd_t ss = xx*trans(yy);
      Map<const Matrix3d> S(&ss[0]);
      JacobiSVD<Matrix3d> sol(S, ComputeFullU|ComputeFullV);
      Matrix3d I = Matrix3d::Identity();
      I(2, 2) = (sol.matrixV()*sol.matrixU().transpose()).determinant();
      Q_[i] = sol.matrixV()*I*sol.matrixU().transpose();
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
  shared_ptr<edge2cell_adjacent> e2c_;
};

class smooth_energy : public surfparam::Functional<double>
{
public:
  typedef zjucad::matrix::matrix<size_t> mati_t;
  typedef zjucad::matrix::matrix<double> matd_t;
  smooth_energy(const mati_t &tris, const matd_t &nods, const double w)
    : tris_(tris), nods_(nods), w_(w) {
    surfparam::cotmatrix(tris_, nods_, 3, &L_);
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
  : max_iter_(5000),
    tolerance_(1e-12),
    lambda_(0.1) {}

frame_field_deform::frame_field_deform(const boost::property_tree::ptree &pt) {
  max_iter_  = pt.get<size_t>("max_iter");
  tolerance_ = pt.get<double>("tolerance");
  lambda_    = pt.get<double>("lambda");
}

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
    matd_t u = nods_(colon(), tris_(1, i))-nods_(colon(), tris_(0, i));
    u /= norm(u);
    matd_t v = cross(normal(colon(), i), u);
    v /= norm(v);
    std::copy(&u[0], &u[0]+3, &B_(0, 3*i+0));
    std::copy(&v[0], &v[0]+3, &B_(0, 3*i+1));
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
//    F_.col(2*fid+0) = (int)(F_.col(2*fid+0).norm()+0.5)*B_.col(3*fid+0);
//    F_.col(2*fid+1) = (int)(F_.col(2*fid+1).norm()+0.5)*B_.col(3*fid+1);
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
    for (size_t k = 0; k < 2; ++k) {
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
  surfparam::rc_vector_row(dx, g2l_, DX);
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
    e_ = std::make_shared<surfparam::energy_t<double>>(buff_);
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
  Map<VectorXd> X(&_nods_[0], _nods_.size());
  for (size_t iter = 0; iter < max_iter_; ++iter) {
    cout << "[info] iteration " << iter << endl;
    // query energy value
    double val = 0;
    e_->Val(&X[0], &val);
    cout << "[info] energy: " << val << "\n";

    // assemble rhs
    VectorXd rhs = VectorXd::Zero(e_->Nx());
    e_->Gra(&X[0], &rhs[0]);
    rhs = -rhs;
    cout << "[info] gradient norm: " << rhs.norm() << "\n\n";

    VectorXd dx = sol_.solve(rhs);
    ASSERT(sol_.info() == Success);

    double x0norm = X.norm();
    X += dx;
    std::dynamic_pointer_cast<deform_energy>(buff_[DEFORM])
        ->UpdateRotation(&X[0]);

    // convergence test
    if ( dx.norm() <= tolerance_*x0norm ) {
      cout << "[info] converged!\n";
      break;
    }
  }
  return 0;
}

}
