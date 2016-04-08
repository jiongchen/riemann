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

//void diffuse_rotation(const mati_t &tris, const matd_t &vrest, const matd_t &vcurr,
//                      const size_t root_face, const tree_t &g,
//                      vector<Matrix3d> &rot) {
//  if ( rot.size() != tris.size(2) )
//    rot.resize(tris.size(2));
//  matd_t v_rest_root = vrest(colon(), tris(colon(), root_face));
//  matd_t v_curr_root = vcurr(colon(), tris(colon(), root_face));
//  rot[root_face] = compute_tri_rot(&v_rest_root[0], &v_curr_root[0]);

//  queue<size_t> q;
//  unordered_set<size_t> vis;
//  q.push(root_face);
//  while ( !q.empty() ) {
//    const size_t curr_face = q.front();
//    q.pop();
//    if ( vis.find(curr_face) != vis.end() )
//      continue;
//    vis.insert(curr_face);
//    for (size_t e = g.first[curr_face]; e != -1; e = g.next[e]) {
//      const size_t next_face = g.v[e];
//      if ( vis.find(next_face) != vis.end() )
//        continue;
//      mati_t diam;
//      get_edge_diam_elem(curr_face, next_face, tris, diam);
//      matd_t x0 = vrest(colon(), diam), x1 = vcurr(colon(), diam);
//      double theta0 = 0, theta1 = 0;
//      calc_dihedral_angle_(&theta0, &x0[0]);
//      calc_dihedral_angle_(&theta1, &x1[0]);
//      bool need_swap = false;
//      for (size_t k = 0; k < 3; ++k) {
//        if ( diam[1] == tris(k, curr_face) ) {
//          if ( diam[2] == tris((k+1)%3, curr_face) ) {
//            need_swap = true;
//            break;
//          }
//        }
//      }
//      matd_t axis = vcurr(colon(), diam[1])-vcurr(colon(), diam[2]);
//      axis *= need_swap ? -1 : 1;
//      rot[next_face] = axis_angle_rot_mat(&axis[0], theta1-theta0)*rot[curr_face];
//      q.push(next_face);
//    }
//  }
//  ASSERT(vis.size() == tris.size(2));
//}

static void tri2tet(const mati_t &tris, const matd_t &v_tri, mati_t &tets, matd_t &v_tet) {
  tets.resize(4, tris.size(2));
  v_tet.resize(3, v_tri.size(2)+tris.size(2));
  tets(colon(0, 2), colon()) = tris;
  tets(3, colon()) = colon(v_tri.size(2), v_tet.size(2)-1);
  v_tet(colon(), colon(0, v_tri.size(2)-1)) = v_tri;
#pragma omp parallel for
  for (size_t i = 0; i < tris.size(2); ++i) {
    matd_t vert = v_tri(colon(), tris(colon(), i));
    matd_t n = cross(vert(colon(), 1)-vert(colon(), 0), vert(colon(), 2)-vert(colon(), 0));
    v_tet(colon(), i+v_tri.size(2)) = vert(colon(), 0)+n/norm(n);
  }
}

class diffuse_arap_energy : public Functional<double>
{
public:
  diffuse_arap_energy(const mati_t &tris, const matd_t &nods, const double w=1.0)
    : w_(w) {
    tri2tet(tris, nods, tets_, nods_);
    dim_ = nods_.size();
    vol_.resize(tets_.size(2));
    D_.resize(tets_.size(2));
    R_.resize(tets_.size(2));
#pragma omp parallel for
    for (size_t i = 0; i < tets_.size(2); ++i) {
      matd_t basis = nods_(colon(), tets_(colon(1, 3), i))-nods_(colon(), tets_(0, i))*ones<double>(1, 3);
      Map<Matrix3d> ba(&basis[0]);
      vol_[i] = fabs(ba.determinant())/6.0;
      D_[i] = ba.inverse();
      R_[i] = Matrix3d::Identity();
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
  void SetRotation(const vector<Matrix3d> &R) {
#pragma omp parallel for
    for (size_t i = 0; i < R.size(); ++i)
      R_[i] = R[i];
  }
private:
  mati_t tets_;
  matd_t nods_;
  size_t dim_;
  const double w_;
  vector<double> vol_;
  vector<Matrix3d> R_;
  vector<Matrix3d> D_;
};
//==============================================================================
void diffuse_arap_encoder::calc_delta_angle(const mati_t &tris, const matd_t &prev, const matd_t &curr,
                                            const tree_t &g, const size_t root_face, vector<double> &da) {
  queue<size_t> q;
  unordered_set<size_t> vis;
  q.push(root_face);
  while ( !q.empty() ) {
    const size_t curr_face = q.front();
    q.pop();
    if ( vis.find(curr_face) != vis.end() ) {
      cout << "[Error] not a tree\n";
      ASSERT(0);
    }
    vis.insert(curr_face);
    for (size_t e = g.first[curr_face]; e != -1; e = g.next[e]) {
      const size_t next_face = g.v[e];
      if ( vis.find(next_face) != vis.end() )
        continue;
      mati_t diam;
      get_edge_diam_elem(curr_face, next_face, tris, diam);
      matd_t x0 = prev(colon(), diam), x1 = curr(colon(), diam);
      double theta0 = 0, theta1 = 0;
      calc_dihedral_angle_(&theta0, &x0[0]);
      calc_dihedral_angle_(&theta1, &x1[0]);
      da.push_back(theta1-theta0);
      q.push(next_face);
    }
  }
  ASSERT(da.size() == tris.size(2)-1);
}

diffuse_arap_decoder::diffuse_arap_decoder(const mati_t &tris, const matd_t &nods)
  : tris_(tris), nods_(nods) {
  energy_ = make_shared<diffuse_arap_energy>(tris, nods);
  dim_ = energy_->Nx();
  X_ = VectorXd::Zero(dim_);
  R_.resize(tris.size(2));
}

int diffuse_arap_decoder::estimate_rotation(const matd_t &prev, const tree_t &g,
                                            const size_t root_face, const vector<double> &da) {
  matd_t root_rest = nods_(colon(), tris_(colon(), root_face));
  matd_t root_prev = prev(colon(), tris_(colon(), root_face));
  R_[root_face] = compute_tri_rot(&root_rest[0], &root_prev[0]);

  queue<size_t> q;
  unordered_set<size_t> vis;
  q.push(root_face);
  size_t cnt = 0;
  while ( !q.empty() ) {
    const size_t curr_face = q.front();
    q.pop();
    if ( vis.find(curr_face) != vis.end() ) {
      cerr << "[Error] not a tree\n";
      ASSERT(0);
    }
    vis.insert(curr_face);
    for (size_t e = g.first[curr_face]; e != -1; e = g.next[e]) {
      const size_t next_face = g.v[e];
      if ( vis.find(next_face) != vis.end() )
        continue;
      mati_t diam;
      get_edge_diam_elem(curr_face, next_face, tris_, diam);
      bool need_swap = false;
      for (size_t k = 0; k < 3; ++k) {
        if ( diam[1] == tris_(k, curr_face) ) {
          if ( diam[2] == tris_((k+1)%3, curr_face) ) {
            need_swap = true;
            break;
          }
        }
      }
      matd_t axis = itr_matrix<const double *>(3, 3, R_[curr_face].data())*(prev(colon(), diam[1])-prev(colon(), diam[2]));
      axis *= need_swap ? -1 : 1;
      R_[next_face] = axis_angle_rot_mat(&axis[0], da[cnt++])*R_[curr_face];
      q.push(next_face);
    }
  }
  ASSERT(vis.size() == tris_.size(2));
  energy_->SetRotation(R_);
  return 0;
}

int diffuse_arap_decoder::pin_down_vert(const size_t id, const double *pos) {
  fixDoF_.insert(3*id+0);
  fixDoF_.insert(3*id+1);
  fixDoF_.insert(3*id+2);
  std::copy(pos, pos+3, X_.data()+3*id);
}

int diffuse_arap_decoder::solve(matd_t &curr) {
  VectorXd g = VectorXd::Zero(dim_); {
    energy_->Gra(&X_[0], g.data());
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
  X_ += DX;
  std::copy(X_.data(), X_.data()+nods_.size(), &curr[0]);
  return 0;
}

}
