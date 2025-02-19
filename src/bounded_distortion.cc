#include "bounded_distortion.h"

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <unordered_map>

#include "config.h"
#include "timer.h"

using namespace std;
using namespace Eigen;
using namespace zjucad::matrix;

namespace riemann {

static SimplicialCholesky<SparseMatrix<double>> ldlt_solver;

class bd_pos_constraint
{
public:
  bd_pos_constraint(const matd_t &nods)
    : dim_(nods.size()) {}
  size_t nx() const { return dim_; }
  size_t nf() const { return 3*fixed_.size(); }
  int rhs(double *val) const {
    RETURN_WITH_COND_TRUE(fixed_.empty());
    Map<VectorXd> fx(val, nf());
    size_t i = 0;
    for (auto &elem : fixed_) {
      fx.segment<3>(3*i++) = elem.second;
    }
    return 0;
  }
  int lhs(const size_t off, vector<Triplet<double>> *jac) const {
    RETURN_WITH_COND_TRUE(fixed_.empty());
    size_t i = 0;
    for (auto &elem : fixed_) {
      jac->push_back(Triplet<double>(off+i++, 3*elem.first+0, 1.0));
      jac->push_back(Triplet<double>(off+i++, 3*elem.first+1, 1.0));
      jac->push_back(Triplet<double>(off+i++, 3*elem.first+2, 1.0));
    }
    return 0;
  }
  int pin(const size_t id, const double *pos) {
    if ( id < 0 || id >= nx()/3 )
      return __LINE__;
    fixed_[id] = Vector3d(pos);
    return 0;
  }
private:
  const size_t dim_;
  unordered_map<size_t, Vector3d> fixed_;
};

static inline void AddDiagBlock3d(const size_t i, const size_t j, const double val, vector<Triplet<double>> &trips) {
  trips.push_back(Triplet<double>(3*i+0, 3*j+0, val));
  trips.push_back(Triplet<double>(3*i+1, 3*j+1, val));
  trips.push_back(Triplet<double>(3*i+2, 3*j+2, val));
}

int calc_tet_base_inv(const mati_t &tets, const matd_t &nods, matd_t &binv) {
  binv.resize(9, tets.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < tets.size(2); ++i) {
    matd_t base = nods(colon(), tets(colon(1, 3), i))-nods(colon(), tets(0, i))*ones<double>(1, 3);
    Map<Matrix3d>(&binv(0, i)) = Map<Matrix3d>(&base[0]).inverse();
  }
  return 0;
}

int calc_tet_df_map(const mati_t &tets, const matd_t &binv, Eigen::SparseMatrix<double> *T) {
  vector<Triplet<double>> trips;
  for (size_t i = 0; i < tets.size(2); ++i) {
    AddDiagBlock3d(3*i+0, tets(0, i), -sum(binv(colon(0, 2), i)), trips);
    AddDiagBlock3d(3*i+0, tets(1, i), binv(0, i), trips);
    AddDiagBlock3d(3*i+0, tets(2, i), binv(1, i), trips);
    AddDiagBlock3d(3*i+0, tets(3, i), binv(2, i), trips);

    AddDiagBlock3d(3*i+1, tets(0, i), -sum(binv(colon(3, 5), i)), trips);
    AddDiagBlock3d(3*i+1, tets(1, i), binv(3, i), trips);
    AddDiagBlock3d(3*i+1, tets(2, i), binv(4, i), trips);
    AddDiagBlock3d(3*i+1, tets(3, i), binv(5, i), trips);

    AddDiagBlock3d(3*i+2, tets(0, i), -sum(binv(colon(6, 8), i)), trips);
    AddDiagBlock3d(3*i+2, tets(1, i), binv(6, i), trips);
    AddDiagBlock3d(3*i+2, tets(2, i), binv(7, i), trips);
    AddDiagBlock3d(3*i+2, tets(3, i), binv(8, i), trips);
  }
  T->resize(9*tets.size(2), 3*(max(tets)+1));
  T->reserve(trips.size());
  T->setFromTriplets(trips.begin(), trips.end());
  return 0;
}

bd_solver::bd_solver(const mati_t &tets, const matd_t &nods, const bd_args &args)
  : tets_(tets), nods_(nods), dim_(nods.size()), lift_dim_(9*tets_.size(2)), args_(args) {
  matd_t binv;
  calc_tet_base_inv(tets_, nods_, binv);
  calc_tet_df_map(tets_, binv, &T_);
  linc_ = make_shared<bd_pos_constraint>(nods_);
}

void bd_solver::set_bound(const double K) {
  args_.K = K;
}

int bd_solver::pin_down_vert(const size_t id, const double *pos) {
  return linc_->pin(id, pos);
}

int bd_solver::prefactorize() {
  cout << "[info] prefactorization......";
  ldlt_solver.setMode(SimplicialCholeskyLDLT);
  SparseMatrix<double> M(dim_+linc_->nf(), dim_+linc_->nf()); {
    vector<Triplet<double>> trips;
    SparseMatrix<double> TtT = T_.transpose()*T_;
    for (size_t j = 0; j < TtT.outerSize(); ++j) {
      for (SparseMatrix<double>::InnerIterator it(TtT, j); it; ++it) {
        trips.push_back(Triplet<double>(it.row(), it.col(), it.value()));
      }
    } {
      vector<Triplet<double>> temp;
      linc_->lhs(dim_, &temp);
      for (auto &elem : temp) {
        trips.push_back(Triplet<double>(elem.row(), elem.col(), elem.value()));
        trips.push_back(Triplet<double>(elem.col(), elem.row(), elem.value()));
      }
    }
    M.reserve(trips.size());
    M.setFromTriplets(trips.begin(), trips.end());
  }
  ldlt_solver.compute(M);
  ASSERT(ldlt_solver.info() == Success);
  cout << "...done\n";
  return 0;
}

int bd_solver::optimize(double *init_x) const {
  int rtn = 0;
  switch ( args_.method ) {
    case 0: rtn = solve(init_x); break;
    case 1: rtn = alter_solve(init_x); break;
    case 2: rtn = alter_solve_chebyshev(init_x); break;
    default: rtn = __LINE__; break;
  }
  return rtn;
}

int bd_solver::solve(double *initX) const {
  cout << "[info] solve\n";
  Map<VectorXd> X(initX, dim_);
  const size_t cdim = linc_->nf();
  VectorXd z(lift_dim_), Pz(lift_dim_), n(lift_dim_), eta = VectorXd::Zero(dim_+cdim);
  VectorXd u(dim_+cdim+1), rhs(dim_+cdim+1), temp0, temp1;
  rhs.setZero();
  linc_->rhs(&rhs[dim_]);
  VectorXd post_step(lift_dim_);
  // linearly constrained least squares
  high_resolution_timer clk;
  clk.start();
  for (size_t iter = 0; iter < args_.maxiter; ++iter) {
    z = T_*X;
    euclidean_proj(&z[0], &Pz[0]);
    n = z-Pz;
    if ( iter % 1 == 0 ) {
      cout << "\t@iter " << iter << " normal norm: " << n.norm() << endl;
      clk.stop();
      clk.log();
    }
    if ( n.norm() < args_.tolerance ) {
      cout << "\t@converged after " << iter << " iterations\n";
      break;
    }
    // linear solve
    eta.head(dim_) = T_.transpose()*n;
    rhs.head(dim_) = T_.transpose()*T_*X;
    rhs[rhs.size()-1] = n.dot(Pz);
    temp0 = ldlt_solver.solve(rhs.head(dim_+cdim));
    ASSERT(ldlt_solver.info() == Success);
    temp1 = ldlt_solver.solve(eta);
    ASSERT(ldlt_solver.info() == Success);
    u[u.size()-1] = (-rhs[rhs.size()-1]+eta.dot(temp0))/eta.dot(temp1);
    u.head(dim_+cdim) = temp0-u[u.size()-1]*temp1;
    X = u.head(dim_);
    post_step = T_*X-Pz;
    if ( iter % 1 == 0 ) {
      cout << "\t@prev step size: " << n.norm() << endl;
      cout << "\t@post step size: " << post_step.norm() << endl;
      cout << "\t@turning angle: " << acos(post_step.dot(n)/(n.norm()*post_step.norm()))/M_PI*180 << "\n\n";
    }
  }
  return 0;
}

int bd_solver::alter_solve(double *initX) const {
  cout << "[info] alternating solve\n";
  Map<VectorXd> X(initX, dim_);
  const size_t cdim = linc_->nf();
  VectorXd z(lift_dim_), Pz(lift_dim_), n(lift_dim_), u(dim_+cdim), rhs(dim_+cdim);
  linc_->rhs(&rhs[dim_]);
  VectorXd post_step(lift_dim_);
  // solve KKT
  double prev_err_norm = 1;
  high_resolution_timer clk;
  clk.start();
  for (size_t iter = 0; iter < args_.maxiter; ++iter) {
    z = T_*X;
    euclidean_proj(&z[0], &Pz[0]);
    n = z-Pz;
    double curr_err_norm = n.norm();
    if ( curr_err_norm < args_.tolerance ) {
      cout << "\t@CONVERGED after " << iter << " iterations\n";
      break;
    }
    if ( iter % 1 == 0 ) {
      cout << "\t@iter " << iter << " normal norm: " << n.norm() << endl;
      cout << "\t@spectral radius: " << curr_err_norm/prev_err_norm << endl;
      clk.stop();
//      clk.log();
    }
    prev_err_norm = curr_err_norm;
    rhs.head(dim_) = T_.transpose()*Pz;
    u = ldlt_solver.solve(rhs);
    ASSERT(ldlt_solver.info() == Success);
    X = u.head(dim_);
    post_step = T_*X-Pz;
    if ( iter % 1 == 0 ) {
      cout << "\t@prev step size: " << n.norm() << endl;
      cout << "\t@post step size: " << post_step.norm() << endl;
      cout << "\t@turning angle: " << acos(post_step.dot(n)/(n.norm()*post_step.norm()))/M_PI*180 << "\n\n";
    }
  }
  return 0;
}

int bd_solver::alter_solve_chebyshev(double *initX) const {
  // x_{k+1} = \omega(\gamma(\hat x_{k+1}-x_k)+x_k-x_{k-1})+x_{k-1}
  cout << "[INFO] alternating solve using chebyshev\n";
  Map<VectorXd> X(initX, dim_);
  const size_t cdim = linc_->nf();
  VectorXd curr_x(dim_), prev_x = X, z(lift_dim_), Pz(lift_dim_), n(lift_dim_), u(dim_+cdim), rhs(dim_+cdim);
  linc_->rhs(&rhs[dim_]);

  static const size_t S = 15;
  static const double rho = args_.sr, gamma = 0.75;
  double omega;
  for (size_t iter = 0; iter < args_.maxiter; ++iter) {
    z = T_*X;
    euclidean_proj(&z[0], &Pz[0]);
    n = z-Pz;
    if ( iter % 1 == 0 ) {
      cout << "\t@iter " << iter << " normal norm: " << n.norm() << endl;
    }
    if ( n.norm() < args_.tolerance ) {
      cout << "\t@CONVERGED after " << iter << " iterations\n";
      break;
    }
    rhs.head(dim_) = T_.transpose()*Pz;
    u = ldlt_solver.solve(rhs);
    ASSERT(ldlt_solver.info() == Success);
    curr_x = X;
    X = u.head(dim_);
    if ( iter < S )
      omega = 1.0;
    else if ( iter == S )
      omega = 2.0/(2.0-rho*rho);
    else
      omega = 4.0/(4.0-rho*rho*omega);
    X = omega*(gamma*(X-curr_x)+curr_x-prev_x)+prev_x;
    prev_x = curr_x;
  }
  return 0;
}

int bd_solver::euclidean_proj(const double *Tx, double *PTx) const {
#pragma omp parallel for
  for (size_t i = 0; i < tets_.size(2); ++i) {
    Map<const Matrix3d> df(Tx+9*i);
    JacobiSVD<Matrix3d> svd(df, ComputeFullU|ComputeFullV);
    Matrix3d U = svd.matrixU(), V = svd.matrixV();
    Vector3d diag = svd.singularValues();
    double detUV = U.determinant()*V.determinant();
    if ( detUV < 0 ) {
      diag[2] *= -1;
      U.col(2) *= -1;
    }
    if ( diag[0] <= args_.K*diag[2] ) {
      Map<Matrix3d>(PTx+9*i) = df;
      continue;
    }
    double t = (args_.K*diag[0]+diag[2])/(1+args_.K*args_.K);
    if ( args_.K*t >= diag[1] && diag[1] >= t ) {
      diag[0] = args_.K*t;
      diag[2] = t;
      Map<Matrix3d>(PTx+9*i) = U*diag.asDiagonal()*V.transpose();
      continue;
    }
    if ( diag[1] > args_.K*t ) {
      double tt = (args_.K*diag[0]+args_.K*diag[1]+diag[2])/(1+2*args_.K*args_.K);
      diag[0] = diag[1] = args_.K*tt;
      diag[2] = tt;
      Map<Matrix3d>(PTx+9*i) = U*diag.asDiagonal()*V.transpose();
      continue;
    }
    if ( t > diag[1] ) {
      double tt = (args_.K*diag[0]+diag[1]+diag[2])/(2+args_.K*args_.K);
      diag[0] = args_.K*tt;
      diag[1] = diag[2] = tt;
      Map<Matrix3d>(PTx+9*i) = U*diag.asDiagonal()*V.transpose();
    }
  }
  return 0;
}

int bd_solver::calc_df_cond_number(const double *x, matd_t &cond) const {
  Map<const VectorXd> X(x, dim_);
  VectorXd z = T_*X;
  cond = zeros<double>(tets_.size(2), 1);
#pragma omp parallel for
  for (size_t i = 0; i < tets_.size(2); ++i) {
    Map<const Matrix3d> df(z.data()+9*i);
    JacobiSVD<Matrix3d> svd(df, ComputeFullU|ComputeFullV);
    cond[i] = svd.singularValues()[0]/svd.singularValues()[2];
  }
  return 0;
}

}
