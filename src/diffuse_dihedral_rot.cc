#include "diffuse_dihedral_rot.h"

#include <queue>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Geometry>

#include "def.h"
#include "dual_graph.h"
#include "config.h"
#include "util.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

namespace riemann {

extern "C" {

void calc_dihedral_angle_(double *val, const double *x);
void tet_arap_(double *val, const double *x, const double *D, const double *R, const double *vol);
void tet_arap_jac_(double *jac, const double *x, const double *D, const double *R, const double *vol);
void tet_arap_hes_(double *hes, const double *x, const double *D, const double *R, const double *vol);

}

static Matrix3d compute_tri_rot(const double *rest, const double *defo) {
  Map<const Matrix3d> X(rest), Y(defo);
  Matrix3d Dm, Ds;
  Dm.col(0) = X.col(1)-X.col(0); Dm.col(1) = X.col(2)-X.col(0); Dm.col(2) = Dm.col(0).cross(Dm.col(1)).normalized();
  Ds.col(0) = Y.col(1)-Y.col(0); Ds.col(1) = Y.col(2)-Y.col(0); Ds.col(2) = Ds.col(0).cross(Ds.col(1)).normalized();
  Matrix3d G = Ds*Dm.inverse();
  JacobiSVD<Matrix3d> svd(G,  ComputeFullU|ComputeFullV);
  return svd.matrixU()*svd.matrixV().transpose();
}

static void get_edge_diam_elem(const size_t f0, const size_t f1, const mati_t &tris, mati_t &diam) {
  // f0 and f1 should be adjacent
  diam.resize(4, 1); {
    mati_t tri0 = tris(colon(), f0), tri1 = tris(colon(), f1);
    vector<size_t> buffer;
    if ( find(tri1.begin(), tri1.end(), tri0[0]) != tri1.end() ) buffer.push_back(tri0[0]);
    if ( find(tri1.begin(), tri1.end(), tri0[1]) != tri1.end() ) buffer.push_back(tri0[1]);
    if ( find(tri1.begin(), tri1.end(), tri0[2]) != tri1.end() ) buffer.push_back(tri0[2]);
    ASSERT(buffer.size() == 2);
    diam[1] = buffer[0]; diam[2] = buffer[1];
  }
  bool need_swap = true;
  for (size_t k = 0; k < 3; ++k) {
    if ( diam[1] == tris(k, f0) ) {
      if ( diam[2] != tris((k+1)%3, f0) )
        need_swap = false;
    }
  }
  if ( need_swap )
    swap(diam[1], diam[2]);
  diam[0] = sum(tris(colon(), f0))-sum(diam(colon(1, 2)));
  diam[3] = sum(tris(colon(), f1))-sum(diam(colon(1, 2)));
}

static inline Matrix3d axis_angle_rot_mat(const double *axis, const double angle) {
  return AngleAxisd(angle, Vector3d(axis).normalized()).toRotationMatrix();
}

void diffuse_rotation(const mati_t &tris, const matd_t &vrest, const matd_t &vcurr,
                      const size_t root_face, const tree_t &g,
                      vector<Matrix3d> &rot) {
  if ( rot.size() != tris.size(2) )
    rot.resize(tris.size(2));
  matd_t v_rest_root = vrest(colon(), tris(colon(), root_face));
  matd_t v_curr_root = vcurr(colon(), tris(colon(), root_face));
  rot[root_face] = compute_tri_rot(&v_rest_root[0], &v_curr_root[0]);

  queue<size_t> q;
  unordered_set<size_t> vis;
  q.push(root_face);
  while ( !q.empty() ) {
    const size_t curr_face = q.front();
    q.pop();
    if ( vis.find(curr_face) != vis.end() )
      continue;
    vis.insert(curr_face);
    for (size_t e = g.first[curr_face]; e != -1; e = g.next[e]) {
      const size_t next_face = g.v[e];
      if ( vis.find(next_face) != vis.end() )
        continue;
      mati_t diam;
      get_edge_diam_elem(curr_face, next_face, tris, diam);
      matd_t x0 = vrest(colon(), diam), x1 = vcurr(colon(), diam);
      double theta0 = 0, theta1 = 0;
      calc_dihedral_angle_(&theta0, &x0[0]);
      calc_dihedral_angle_(&theta1, &x1[0]);
      matd_t axis = vcurr(colon(), diam[1])-vcurr(colon(), diam[2]);
      rot[next_face] = axis_angle_rot_mat(&axis[0], theta1-theta0)*rot[curr_face];
      q.push(next_face);
    }
  }
  ASSERT(vis.size() == tris.size(2));
}

class diffuse_arap_energy : public Functional<double>
{
public:
  diffuse_arap_energy(const mati_t &tets, const matd_t &nods, const vector<Matrix3d> &R, const double w=1.0)
    : tets_(tets), nods_(nods), dim_(nods.size()), R_(R), w_(w) {
    vol_.resize(tets_.size(2));
    D_.resize(tets_.size(2));
#pragma omp parallel for
    for (size_t i = 0; i < tets_.size(2); ++i) {
      matd_t basis = nods_(colon(), tets_(colon(1, 3), i))-nods_(colon(), tets_(0, i))*ones<double>(1, 3);
      Map<Matrix3d> ba(&basis[0]);
      vol_[i] = fabs(ba.determinant())/6.0;
      D_[i] = ba.inverse();
    }
  }
  size_t Nx() const {
    return nods_.size();
  }
  int Val(const double *x, double *val) const {
    itr_matrix<const double *> X(3, dim_/3, x);
    for (size_t i = 0; i < tets_.size(2); ++i) {
      matd_t vert = X(colon(), tets_(colon(), i));
      double value = 0;
      tet_arap_(&value, &vert[0], D_[i].data(), R_[i].data(), &vol_[i]);
      *val += w_*value;
    }
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    itr_matrix<const double *> X(3, dim_/3, x);
    itr_matrix<double *> G(3, dim_/3, gra);
    for (size_t i = 0; i < tets_.size(2); ++i) {
      matd_t vert = X(colon(), tets_(colon(), i));
      matd_t grad = zeros<double>(3, 4);
      tet_arap_jac_(&grad[0], &vert[0], D_[i].data(), R_[i].data(), &vol_[i]);
      G(colon(), tets_(colon(), i)) += w_*grad;
    }
    return 0;
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    for (size_t i = 0; i < tets_.size(2); ++i) {
      matd_t H = zeros<double>(12, 12);
      tet_arap_hes_(&H[0], nullptr, D_[i].data(), R_[i].data(), &vol_[i]);
      for (size_t p = 0; p < 12; ++p) {
        for (size_t q = 0; q < 12; ++q) {
          const size_t I = 3*tets_(p/3, i)+p%3;
          const size_t J = 3*tets_(q/3, i)+q%3;
          hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
        }
      }
    }
    return 0;
  }
private:
  const mati_t &tets_;
  const matd_t &nods_;
  const size_t dim_;
  double w_;
  vector<double> vol_;
  const vector<Matrix3d> &R_;
  vector<Matrix3d> D_;
};
//==============================================================================
diffuse_arap_solver::diffuse_arap_solver(const mati_t &tets, const matd_t &nods, const vector<Matrix3d> &R)
  : dim_(nods.size()) {
  energy_ = make_shared<diffuse_arap_energy>(tets, nods, R);
}

int diffuse_arap_solver::pin_down_vert(const size_t id, const double *pos, double *x) {
  fixDoF_.insert(3*id+0);
  fixDoF_.insert(3*id+1);
  fixDoF_.insert(3*id+2);
  std::copy(pos, pos+3, x+3*id);
}

int diffuse_arap_solver::solve(double *x) const {
  Map<VectorXd> X(x, dim_);
  VectorXd g = VectorXd::Zero(dim_); {
    energy_->Gra(&x[0], g.data());
    g *= -1;
  }
  SparseMatrix<double> H(dim_, dim_); {
    vector<Triplet<double>> trips;
    energy_->Hes(nullptr, &trips);
    H.setFromTriplets(trips.begin(), trips.end());
  }
  vector<size_t> g2l(dim_);
  size_t cnt = 0;
  for (size_t i = 0; i < dim_; ++i) {
    if ( fixDoF_.find(i) != fixDoF_.end() )
      g2l[i] = -1;
    else
      g2l[i] = cnt++;
  }
  SimplicialCholesky<SparseMatrix<double>> solver;
  rm_spmat_col_row(H, g2l);
  rm_vector_row(g, g2l);
  solver.compute(H);
  ASSERT(solver.info() == Success);
  VectorXd dx = solver.solve(g);
  ASSERT(solver.info() == Success);
  VectorXd DX = VectorXd::Zero(dim_);
  rc_vector_row(dx, g2l, DX);
  X += DX;
  return 0;
}

}
