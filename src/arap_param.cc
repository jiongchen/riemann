#include "arap_param.h"

#include <iostream>
#include <zjucad/matrix/itr_matrix.h>
#include <Eigen/SVD>

#include "def.h"
#include "geometry_extend.h"
#include "util.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

namespace riemann {

static void calc_tris_cot_value(const mati_t &tris, const matd_t &nods, matd_t &cotv) {
  cotv.resize(3, tris.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < tris.size(2); ++i) {
    matd_t vert = nods(colon(), tris(colon(), i));
    cotv(0, i) = cal_cot_val(&vert(0, 2), &vert(0, 0), &vert(0, 1));
    cotv(1, i) = cal_cot_val(&vert(0, 0), &vert(0, 1), &vert(0, 2));
    cotv(2, i) = cal_cot_val(&vert(0, 1), &vert(0, 2), &vert(0, 0));
  }
}

class arap_param_energy : public Functional<double>
{
public:
  arap_param_energy(const mati_t &tris, const matd_t &nods, const double w=1.0)
    : w_(w), tris_(tris), dim_(2*nods.size(2)), R_(tris.size(2)){
    matd_t origin, axis;
    calc_face_local_frame(tris, nods, origin, axis);
    calc_local_uv(tris, nods, origin, axis, luv_);
    calc_tris_cot_value(tris, nods, cotv_);
    for_each(R_.begin(), R_.end(), [](Matrix2d &t){ t.setIdentity(); });
  }
  size_t Nx() const {
    return dim_;
  }
  int Val(const double *x, double *val) const {
    RETURN_WITH_COND_TRUE(w_ == 0.0);
    itr_matrix<const double *> X(2, dim_/2, x);
    for (size_t i = 0; i < tris_.size(2); ++i) {
      itr_matrix<const double *> R(2, 2, R_[i].data());
      for (size_t j = 0; j < 3; ++j) {
        const size_t p = (j+1)%3, q = (j+2)%3, gp = tris_(p, i), gq = tris_(q, i);
        matd_t d = X(colon(), gp)-X(colon(), gq)-R*(luv_(colon(2*p, 2*p+1), i)-luv_(colon(2*q, 2*q+1), i));
        *val += 0.5*w_*cotv_(j, i)*dot(d, d);
      }
    }
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    RETURN_WITH_COND_TRUE(w_ == 0.0);
    itr_matrix<const double *> X(2, dim_/2, x);
    itr_matrix<double *> G(2, dim_/2, gra);
    for (size_t i = 0; i < tris_.size(2); ++i) {
      itr_matrix<const double *> R(2, 2, R_[i].data());
      for (size_t j = 0; j < 3; ++j) {
        const size_t p = (j+1)%3, q = (j+2)%3, gp = tris_(p, i), gq = tris_(q, i);
        matd_t d = X(colon(), gp)-X(colon(), gq)-R*(luv_(colon(2*p, 2*p+1), i)-luv_(colon(2*q, 2*q+1), i));
        G(colon(), gp) += w_*cotv_(j, i)*d;
        G(colon(), gq) -= w_*cotv_(j, i)*d;
      }
    }
    return 0;
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    RETURN_WITH_COND_TRUE(w_ == 0.0);
    for (size_t i = 0; i < tris_.size(2); ++i) {
      for (size_t j = 0; j < 3; ++j) {
        const size_t p = (j+1)%3, q = (j+2)%3, gp = tris_(p, i), gq = tris_(q, i);
        const double v = w_*cotv_(j, i);
        add_diag_block<double, 2>(gp, gp, v, hes);
        add_diag_block<double, 2>(gp, gq, -v, hes);
        add_diag_block<double, 2>(gq, gp, -v, hes);
        add_diag_block<double, 2>(gq, gq, v ,hes);
      }
    }
    return 0;
  }
  void LocalSolve(const double *x) {
    Map<const MatrixXd> X(x, 2, dim_/2);
    Map<const MatrixXd> U(&luv_[0], luv_.size(1), luv_.size(2));
#pragma omp parallel for
    for (size_t i = 0; i < tris_.size(2); ++i) {
      Matrix2d S = Matrix2d::Zero();
      for (size_t j = 0; j < 3; ++j) {
        const size_t p = (j+1)%3, q = (j+2)%3, gp = tris_(p, i), gq = tris_(q, i);
        S += cotv_(j, i)*(X.col(gp)-X.col(gq))*(U.block<2, 1>(2*p, i)-U.block<2, 1>(2*q, i)).transpose();
      }
      JacobiSVD<Matrix2d> svd(S, ComputeFullU|ComputeFullV);
      R_[i] = svd.matrixU()*svd.matrixV().transpose();
      if ( R_[i].determinant() < 0.0 )
        R_[i].col(1) *= -1;
    }
  }
private:
  const double w_;
  const mati_t &tris_;
  const size_t dim_;
  matd_t luv_;
  matd_t cotv_;
  vector<Matrix2d> R_;
};

arap_param_solver::arap_param_solver(const mati_t &tris, const matd_t &nods) {
  arap_ = make_shared<arap_param_energy>(tris, nods);
}

int arap_param_solver::precompute() {
  vector<Triplet<double>> trips;
  arap_->Hes(nullptr, &trips);
  LHS_.resize(arap_->Nx(), arap_->Nx());
  LHS_.setFromTriplets(trips.begin(), trips.end());
  solver_.compute(LHS_);
  ASSERT(solver_.info() == Success);
  return 0;
}

int arap_param_solver::solve(double *x0) const {
  Map<VectorXd> x(x0, arap_->Nx());
  VectorXd dx(arap_->Nx());
  for (size_t iter = 0; iter < 2000; ++iter) {
    double value = 0; {
      arap_->Val(&x[0], &value);
      if ( iter % 10 == 0 ) {
        cout << "\t@energy value: " << value << endl;
      }
    }
    // local solve
    dynamic_pointer_cast<arap_param_energy>(arap_)->LocalSolve(&x[0]);
    // global solve
    VectorXd g = VectorXd::Zero(arap_->Nx()); {
      arap_->Gra(&x[0], &g[0]);
      g *= -1;
    }
    // convergence test
    if ( g.norm() <= 1e-8 ) {
      cout << "\t@CONVERGED after " << iter << " iterations\n";
      break;
    }
    dx = solver_.solve(g);
    ASSERT(solver_.info() == Success);
    x += dx;
  }
  return 0;
}

}
