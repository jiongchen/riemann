#include "lscm.h"

#include <jtflib/mesh/mesh.h>

#include "energy.h"
#include "util.h"

using namespace std;
using namespace Eigen;
using namespace zjucad::matrix;

namespace surfparam {

lscm_param::lscm_param(const mati_t &tris, const matd_t &nods)
  : tris_(tris) {
  buff_.push_back(shared_ptr<Functional<double>>(new dirichlet_energy(tris_, nods)));
  buff_.push_back(shared_ptr<Functional<double>>(new param_area(tris_, nods)));
  try {
    conformal_energy_.reset(new energy_t<double>(buff_));
  } catch ( exception &e ) {
    cerr << "# exception: " << e.what() << endl;
    exit(0);
  }
  uv_ = VectorXd::Zero(2*nods.size(2));
  get_boundary_loop(tris_, bnd_);
}

int lscm_param::get_boundary_loop(const mati_t &cell, vector<size_t> &bnd) {
  return 0;
}

void lscm_param::set_fixed_bnd_vert(const size_t id, const double *x) {
  fixed_dofs_.insert(2*id+0);
  fixed_dofs_.insert(2*id+1);
  uv_[2*id+0] = x[0];
  uv_[2*id+1] = x[1];
}

int lscm_param::apply() {
  g2l_.resize(uv_.rows());
  size_t ptr = 0;
  for (size_t i = 0; i < g2l_.size(); ++i) {
    if ( fixed_dofs_.find(i) != fixed_dofs_.end() )
      g2l_[i] = -1;
    else
      g2l_[i] = ptr++;
  }

  const size_t dim = conformal_energy_->Nx();
  ///> assemble LHS
  vector<Triplet<double>> trips;
  conformal_energy_->Hes(uv_.data(), &trips);
  SparseMatrix<double> K(dim, dim);
  K.reserve(trips.size());
  K.setFromTriplets(trips.begin(), trips.end());

  ///> assemble rhs
  VectorXd rhs = VectorXd::Zero(dim);
  conformal_energy_->Gra(uv_.data(), rhs.data());
  rhs = -rhs;

  if ( !fixed_dofs_.empty() ) {
    rm_spmat_col_row(K, g2l_);
    rm_vector_row(rhs, g2l_);
  }

  SimplicialCholesky<SparseMatrix<double>> sol;
  sol.compute(K);
  if ( sol.info() != Success ) {
    cerr << "# info: failed to precompute system matrix\n";
    return __LINE__;
  }
  VectorXd delta = sol.solve(rhs);
  if ( sol.info() != Success ) {
    cerr << "# info: failed to solve for rhs\n";
    return __LINE__;
  }

  VectorXd dx = VectorXd::Zero(dim);
  if ( !fixed_dofs_.empty() ) {
    rc_vector_row(delta, g2l_, dx);
  } else {
    dx = delta;
  }
  uv_ += dx;

  double value = 0;
  conformal_energy_->Val(uv_.data(), &value);
  cout << "# info: energy: " << value << endl;

  return 0;
}

int lscm_param::apply_spetral() {
  const size_t dim = conformal_energy_->Nx();
  ///> assemble Lc
  vector<Triplet<double>> trips;
  conformal_energy_->Hes(uv_.data(), &trips);
  SparseMatrix<double> K(dim, dim);
  K.reserve(trips.size());
  K.setFromTriplets(trips.begin(), trips.end());

  ///> construct e
  MatrixXd e = MatrixXd::Zero(dim, 2);
  for (size_t i = 0; i < dim/2; ++i) {
    e(2*i+0, 0) = 1.0;
    e(2*i+1, 1) = 1.0;
  }

  ///> assemble B
  vector<Triplet<double>> btrips;

  ///> assemble B-1/Vb*(Be)(Be)^T

  ///> solve generalized eigenvalue problem
  ///> (B-1/Vb(Be)(Be)^T)*u = miu*Lc*u
  return 0;
}

int lscm_param::get_param_mesh(mati_t *param_tris, matd_t *param_nods) {
  if ( param_tris == nullptr || param_nods == nullptr )
    return __LINE__;
  *param_tris = tris_;
  *param_nods = zeros<double>(3, uv_.size()/2);
#pragma omp parallel for
  for (size_t i = 0; i < param_nods->size(2); ++i) {
    (*param_nods)(0, i) = uv_[2*i+0];
    (*param_nods)(2, i) = uv_[2*i+1];
  }
  return 0;
}

}
