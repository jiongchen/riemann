#include "arap_deform.h"

#include <zjucad/matrix/matrix.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>

using namespace std;
using namespace Eigen;
using namespace zjucad::matrix;
using namespace jtf::mesh;

namespace core {

template <typename T, size_t dim = 3>
void add_diag_block(const size_t row, const size_t col, const T val, vector<Triplet<T>> *mat) {
  const size_t row_offset = dim * row;
  const size_t col_offset = dim * col;
  for (size_t i = 0; i < dim; ++i)
    mat->push_back(Triplet<T>(row_offset+i, col_offset+i, val));
}

template <typename T>
void rm_spmat_col_row(SparseMatrix<T> &A, const vector<size_t> &g2l) {
  size_t new_size = 0;
  for (size_t i = 0; i < g2l.size(); ++i) {
    if ( g2l[i] != -1)
      ++new_size;
  }
  std::vector<Eigen::Triplet<T>> trips;
  for (size_t j = 0; j < A.outerSize(); ++j) {
    for (typename Eigen::SparseMatrix<T>::InnerIterator it(A, j); it; ++it) {
      if ( g2l[it.row()] != -1 && g2l[it.col()] != -1 )
        trips.push_back(Eigen::Triplet<T>(g2l[it.row()], g2l[it.col()], it.value()));
    }
  }
  A.resize(new_size, new_size);
  A.reserve(trips.size());
  A.setFromTriplets(trips.begin(), trips.end());
}

template <typename T>
void rm_spmat_col_row(SparseMatrix<T> &A, const unordered_set<size_t> &idx) {
  vector<size_t> g2l(A.cols());
  size_t ptr = 0;
  for (size_t i = 0; i < g2l.size(); ++i) {
    if ( idx.find(i) != idx.end() )
      g2l[i] = -1;
    else
      g2l[i] = ptr++;
  }
  rm_spmat_col_row<T>(A, g2l);
}

template <typename T>
void rm_vector_row(Matrix<T, -1, 1> &g, const vector<size_t> &g2l) {
  size_t new_size = 0;
  for (size_t i = 0; i < g2l.size(); ++i) {
    if ( g2l[i] != -1 )
      ++new_size;
  }
  Eigen::Matrix<T, -1, 1> sub;
  sub.resize(new_size);
#pragma omp parallel for
  for (size_t i = 0; i < g2l.size(); ++i)
    if ( g2l[i] != -1 )
      sub[g2l[i]] = g[i];
  g = sub;
}

template <typename T>
void up_vector_row(const Matrix<T, -1, 1> &l, const vector<size_t> &g2l, Matrix<T, -1, 1> &g) {
#pragma omp parallel for
  for (size_t i = 0; i < g2l.size(); ++i) {
    if ( g2l[i] != -1 )
      g[i] = l[g2l[i]];
  }
}

template <typename T>
inline T cal_cot_val(const T* a, const T* b, const T* c) {
  Matrix<T, 3, 1> ab(b[0]-a[0], b[1]-a[1], b[2]-a[2]);
  Matrix<T, 3, 1> bc(c[0]-b[0], c[1]-b[1], c[2]-b[2]);
  Matrix<T, 3, 1> ca(a[0]-c[0], a[1]-c[1], a[2]-c[2]);
  return 0.5 * (ab.dot(ab) + bc.dot(bc) - ca.dot(ca)) / ab.cross(bc).norm();
}

class arap_energy
{
public:
  virtual ~arap_energy() {}
  virtual size_t dim() const {}
  virtual int val(const double *x, double *value) const { return __LINE__; }
  virtual int gra(const double *x, double *jac) const { return __LINE__; }
  virtual int hes(const double *x, SparseMatrix<double> *H) const { return __LINE__; }
  virtual int eval_rotation(const double *x) { return __LINE__; }
};

class one_ring_arap_energy : public arap_energy
{
public:
  typedef zjucad::matrix::matrix<size_t> mati_t;
  typedef zjucad::matrix::matrix<double> matd_t;
  one_ring_arap_energy(const mati_t &tris, const matd_t &nods, const double w)
    : tris_(tris), nods_(nods), w_(w) {
    e2c_.reset(edge2cell_adjacent::create(tris_, false));
    wij_ = zeros<double>(e2c_->edges_.size(), 1);
    for (size_t i = 0; i < e2c_->edges_.size(); ++i) {
      size_t pi = e2c_->edges_[i].first;
      size_t qi = e2c_->edges_[i].second;
      pair<size_t, size_t> face = e2c_->query(pi, qi);
      const size_t f[] = {face.first, face.second};
      for (size_t j = 0; j < 2; ++j) {
        if ( f[j] == -1 )
          continue;
        size_t ri = sum(tris_(colon(), f[j]))-pi-qi;
        wij_[i] += 0.5 * cal_cot_val(&nods_(0, pi), &nods_(0, ri), &nods_(0, qi));
      }
    }
    rot_.resize(nods_.size(2));
    shared_ptr<one_ring_point_at_point> p2p(one_ring_point_at_point::create(tris_));
    p2p_.resize(nods_.size(2));
    for (auto &it : p2p->p2p_)
      p2p_[it.first] = it.second;
  }
  size_t dim() const {
    return nods_.size();
  }
  int val(const double *x, double *value) const {
    Map<const MatrixXd> X(x, 3, dim()/3);
    Map<const MatrixXd> X0(&nods_[0], 3, dim()/3);
    for (size_t i = 0; i < e2c_->edges_.size(); ++i) {
      size_t pi = e2c_->edges_[i].first;
      size_t qi = e2c_->edges_[i].second;
      Vector3d eij = X.col(pi)-X.col(qi);
      Vector3d reij = X0.col(pi)-X0.col(qi);
      *value += w_*(wij_[i]*(eij-rot_[pi]*reij).squaredNorm()+
                    wij_[i]*(-eij+rot_[qi]*reij).squaredNorm()
                    );
    }
    return 0;
  }
  int gra(const double *x, double *jac) const {
    Map<const MatrixXd> X(x, 3, dim()/3);
    Map<const MatrixXd> X0(&nods_[0], 3, dim()/3);
    Map<MatrixXd> G(jac, 3, dim()/3);
    size_t i = 0;
    for (auto it : e2c_->edges_) {
      size_t pi = it.first;
      size_t qi = it.second;
      Vector3d eij = X.col(pi)-X.col(qi);
      Vector3d reij = X0.col(pi)-X0.col(qi);
      Vector3d dyi = 2*w_*wij_[i]*(eij-rot_[pi]*reij);
      Vector3d dyj = 2*w_*wij_[i]*(-eij+rot_[qi]*reij);
      G.col(pi) += dyi;
      G.col(qi) -= dyi;
      G.col(qi) += dyj;
      G.col(pi) -= dyj;
      ++i;
    }
    return 0;
  }
  int hes(const double *x, SparseMatrix<double> *H) const {
    vector<Triplet<double>> trips;
    size_t i = 0;
    for (auto &it : e2c_->edges_) {
      size_t pi = it.first;
      size_t qi = it.second;
      double ele = 4*w_*wij_[i];
      add_diag_block(pi, pi, ele, &trips);
      add_diag_block(qi, qi, ele, &trips);
      add_diag_block(pi, qi, -ele, &trips);
      add_diag_block(qi, pi, -ele, &trips);
      ++i;
    }
    H->resize(dim(), dim());
    H->reserve(trips.size());
    H->setFromTriplets(trips.begin(), trips.end());
    return 0;
  }
  int eval_rotation(const double *x) {
    Map<const MatrixXd> X(x, 3, dim()/3);
    Map<const MatrixXd> X0(&nods_[0], 3, dim()/3);
#pragma omp parallel for
    for (size_t pi = 0; pi < p2p_.size(); ++pi) {
      rot_[pi] = Matrix3d::Zero();
      for (auto &qi : p2p_[pi]) {
        size_t eid = e2c_->get_edge_idx(pi, qi);
        Vector3d eij = X.col(pi)-X.col(qi);
        Vector3d reij = X0.col(pi)-X0.col(qi);
        rot_[pi] += wij_[eid]*reij*eij.transpose();
      }
      JacobiSVD<Matrix3d> sol(rot_[pi], ComputeFullU|ComputeFullV);
      rot_[pi] = sol.matrixV()*sol.matrixU().transpose();
    }
    return 0;
  }
private:
  const mati_t tris_;
  const matd_t nods_;
  const double w_;
  matd_t wij_;
  vector<Matrix3d> rot_;
  shared_ptr<edge2cell_adjacent> e2c_;
  vector<vector<size_t>> p2p_;
};

class tri_centric_arap_energy : public arap_energy
{
public:
  typedef zjucad::matrix::matrix<size_t> mati_t;
  typedef zjucad::matrix::matrix<double> matd_t;
  typedef Eigen::Matrix<double, 9, 9> Matrix9d;
  tri_centric_arap_energy(const mati_t &tris, const matd_t &nods, const double w)
    : tris_(tris), nods_(nods), w_(w) {
    /// calculate face normal
    jtf::mesh::cal_face_normal(tris_, nods_, normal_, true);
    /// calculate face area
    area_.resize(1, tris_.size(2));
#pragma omp parallel for
    for (size_t i = 0; i < area_.size(); ++i) {
      matd_t vt = nods_(colon(), tris_(colon(), i));
      area_[i] = jtf::mesh::cal_face_area(vt);
    }
    /// calculate discrete deformation gradient
    /// operator which is piecewise constant
    G_.resize(tris_.size(2));
    //#pragma omp parallel for
    for (size_t i = 0; i < G_.size(); ++i) {
      matd_t Vi = nods_(colon(), tris_(colon(), i));
      G_[i].setZero();
      matd_t Gi(3, 3);
      Gi(colon(), colon(0, 1)) = Vi(colon(), colon(1, 2))
          -Vi(colon(), 0)*ones<double>(1, 2);
      Gi(colon(), 2) = normal_(colon(), i);
      if ( inv(Gi) ) {
        cerr << "#info: inverse fail\n";
      }
      G_[i].block<3, 3>(0, 0) = (-Gi(1, 0)-Gi(0, 0))*Matrix3d::Identity();
      G_[i].block<3, 3>(0, 3) = Gi(0, 0)*Matrix3d::Identity();
      G_[i].block<3, 3>(0, 6) = Gi(1, 0)*Matrix3d::Identity();

      G_[i].block<3, 3>(3, 0) = (-Gi(1, 1)-Gi(0, 1))*Matrix3d::Identity();
      G_[i].block<3, 3>(3, 3) = Gi(0, 1)*Matrix3d::Identity();
      G_[i].block<3, 3>(3, 6) = Gi(1, 1)*Matrix3d::Identity();

      G_[i].block<3, 3>(6, 0) = (-Gi(1, 2)-Gi(0, 2))*Matrix3d::Identity();
      G_[i].block<3, 3>(6, 3) = Gi(0, 2)*Matrix3d::Identity();
      G_[i].block<3, 3>(6, 6) = Gi(1, 2)*Matrix3d::Identity();
    }
    /// allocate space for rotation matrix
    R_.resize(tris_.size(2));
  }
  size_t dim() const {
    return nods_.size();
  }
  int val(const double *x, double *val) const {
    itr_matrix<const double*> X(3, dim()/3, x);
    for (size_t i = 0; i < tris_.size(2); ++i) {
      matd_t vt = X(colon(), tris_(colon(), i));
      Map<const VectorXd> U(&vt[0], 9);
      Map<const VectorXd> R(R_[i].data(), 9);
      *val += w_*area_[i]*(G_[i]*U-R).squaredNorm();
    }
    return 0;
  }
  int gra(const double *x, double *gra) const {
    itr_matrix<const double*> X(3, dim()/3, x);
    itr_matrix<double*> grad(dim(), 1, gra);
    for (size_t i = 0; i < tris_.size(2); ++i) {
      matd_t vt = X(colon(), tris_(colon(), i));
      Map<const VectorXd> U(&vt[0], 9);
      Map<const VectorXd> R(R_[i].data(), 9);
      VectorXd g(9);
      g = 2*w_*area_[i]*(G_[i]*U-R);
      for (size_t j = 0; j < 9; ++j)
        grad[3*tris_(j/3, i)+j%3] += g[j];
    }
    return 0;
  }
  int hes(const double *x, SparseMatrix<double> *H) const {
    vector<Triplet<double>> trip;
    for (size_t i = 0; i < tris_.size(2); ++i) {
      MatrixXd LH = 2*w_*area_[i]*G_[i].transpose()*G_[i];
      for (size_t p = 0; p < 9; ++p) {
        for (size_t q = 0; q < 9; ++q) {
          size_t I = 3*tris_(p/3, i)+p%3;
          size_t J = 3*tris_(q/3, i)+q%3;
          if ( LH(p, q) != 0.0 )
            trip.push_back(Triplet<double>(I, J, LH(p, q)));
        }
      }
    }
    H->resize(dim(), dim());
    H->reserve(trip.size());
    H->setFromTriplets(trip.begin(), trip.end());
    return 0;
  }
  int eval_rotation(const double *x) {
    itr_matrix<const double *> X(3, dim()/3, x);
    matd_t def_normal;
    jtf::mesh::cal_face_normal(tris_, X, def_normal, true);
//#pragma omp parallel for
    for (size_t i = 0; i < tris_.size(2); ++i) {
      matd_t vt = X(colon(), tris_(colon(), i));
      Map<const VectorXd> U(vt.begin(), 9);
      VectorXd def_grad = G_[i]*U;
      Map<Matrix3d> dg(def_grad.data());
//      cout << dg << endl;
//      getchar();
      JacobiSVD<Matrix3d> svd(dg, ComputeFullU|ComputeFullV);
//      cout << sol.singularValues().transpose() << endl;
//      getchar();
      Matrix3d S = svd.matrixU();
      Matrix3d T = svd.matrixV();
      S.col(2) = Vector3d(normal_(0, i), normal_(1, i), normal_(2, i));
      T.col(2) = Vector3d(def_normal(0, i), def_normal(1, i), def_normal(2, i));
      R_[i] = S*T.transpose();
//      static int count = 0;
//      if ( count == 0 && std::fabs(R_[i].squaredNorm()-3) > 1e-12 ) {
//        cout << "id: " << i << endl;
//        cout << R_[i] << endl << endl;
//        getchar();
//      }
    }
    return 0;
  }
private:
  const mati_t tris_;
  const matd_t nods_;
  const double w_;
  matd_t normal_;
  matd_t area_;
  vector<Matrix9d> G_;
  vector<Matrix3d> R_;
};

arap_deform::arap_deform(const mati_t &tris, const matd_t &nods)
  : tris_(tris), nods_(nods) {}

int arap_deform::pre_compute(const vector<size_t> &idx) {
//  e_.reset(new one_ring_arap_energy(tris_, nods_, 1.0));
  e_.reset(new tri_centric_arap_energy(tris_, nods_, 1.0));
  e_->hes(nullptr, &L_);

  fixed_dofs_.clear();
  for (auto &ele : idx) {
    fixed_dofs_.insert(ele*3+0);
    fixed_dofs_.insert(ele*3+1);
    fixed_dofs_.insert(ele*3+2);
  }
  g2l_.resize(L_.cols());
  size_t ptr = 0;
  for (size_t i = 0; i < g2l_.size(); ++i) {
    if ( fixed_dofs_.find(i) != fixed_dofs_.end() )
      g2l_[i] = -1;
    else
      g2l_[i] = ptr++;
  }
  if ( !fixed_dofs_.empty() )
    rm_spmat_col_row(L_, g2l_);
  sol_.compute(L_);
  if ( sol_.info() != Success ) {
    cerr << "[info] prefactorization failed\n";
    return __LINE__;
  }
  return 0;
}

int arap_deform::deformation(double *x) {
  Map<VectorXd> X(x, e_->dim());
  const size_t max_iter = 10000;
  VectorXd Xstar = X;
  VectorXd Dx(e_->dim());
  for (size_t iter = 0; iter < max_iter; ++iter) {
    e_->eval_rotation(Xstar.data());
    if ( iter % 100 == 0 ) {
      double value = 0;
      e_->val(Xstar.data(), &value);
      cout << "[info] iteration " << iter << endl;
      cout << "[info] energy value: " << value << endl << endl;
    }
    VectorXd grad(e_->dim());
    grad.setZero();
    e_->gra(Xstar.data(), grad.data());
    grad = -grad;
    if ( !fixed_dofs_.empty() )
      rm_vector_row(grad, g2l_);
    VectorXd dx = sol_.solve(grad);
    Dx.setZero();
    if ( !fixed_dofs_.empty() )
      up_vector_row(dx, g2l_, Dx);
    else
      Dx = dx;
    double xstar_norm = Xstar.norm();
    Xstar += Dx;
    // convergence test
    if ( Dx.norm() <= 1e-8 * xstar_norm ) {
      printf("# INFO: converged after %zu iteration\n", iter);
      break;
    }
  }
  X = Xstar;
  return 0;
}

}
