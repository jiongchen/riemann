#include "polycube.h"

#include <iostream>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <zjucad/matrix/itr_matrix.h>
#include <jtflib/mesh/mesh.h>
#include <Eigen/CholmodSupport>
#include <zjucad/matrix/io.h>
#include <Eigen/Eigenvalues>

#include "config.h"
#include "def.h"

using namespace std;
using namespace Eigen;
using namespace zjucad::matrix;
using namespace jtf::mesh;

namespace riemann {

extern "C" {

  void surf_normal_align_(double *val, const double *x, const double *eps);
  void surf_normal_align_jac_(double *gra, const double *x, const double *eps);
  void surf_normal_align_hes_(double *hes, const double *x, const double *eps);
  
  void tet_distortion_(double *val, const double *x, const double *D, const double *R);
  void tet_distortion_jac_(double *jac, const double *x, const double *D, const double *R);
  void tet_distortion_hes_(double *hes, const double *x, const double *D, const double *R);

  void triangle_area_(double *val, const double *x);
  void triangle_area_jac_(double *jac, const double *x);
  
}

void area_normal_align_hes(const double *x, const double eps, const int id, double *out);

static inline void calc_def_grad(const double *x, const double *Dm, double *grad) {
  itr_matrix<const double *> X(3, 4, x);
  matd_t Ds = X(colon(), colon(1, 3))-X(colon(), 0)*ones<double>(1, 3);
  itr_matrix<double *>(3, 3, grad) = Ds*itr_matrix<const double *>(3, 3, Dm);
}

class area_preserving_cons : public Functional<double>
{
public:
  area_preserving_cons(const mati_t &surf, const matd_t &nods)
      : surf_(surf), dim_(nods.size()), sum_area_(0) {
    matd_t vert = zeros<double>(3, 3); double area = 0;
    for (size_t i = 0; i < surf_.size(2); ++i) {
      vert = nods(colon(), surf(colon(), i));
      triangle_area_(&area, &vert[0]);
      sum_area_ += area;
    }
    cout << "[Info] surface area: " << sum_area_ << endl;
  }
  size_t Nx() const {
    return dim_;
  }
  int Val(const double *x, double *val) const {
    itr_matrix<const double *> X(3, dim_/3, x);
    matd_t vert = zeros<double>(3, 3); double value = 0;
    double total = 0;
    for (size_t i = 0; i < surf_.size(2); ++i) {
      vert = X(colon(), surf_(colon(), i));
      triangle_area_(&value, &vert[0]);
      total += value;
    }
    *val += total-sum_area_;
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    itr_matrix<const double *> X(3, dim_/3, x);
    itr_matrix<double *> G(3, dim_/3, gra);
    matd_t vert = zeros<double>(3, 3), g = zeros<double>(3, 3);
    for (size_t i = 0; i < surf_.size(2); ++i) {
      vert = X(colon(), surf_(colon(), i));
      triangle_area_jac_(&g[0], &vert[0]);
      G(colon(), surf_(colon(), i)) += g;
    }
    return 0;
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    return __LINE__;
  }
private:
  const mati_t &surf_;
  const size_t dim_;
  double sum_area_;
};

class surf_normal_align_energy : public Functional<double>
{
public:
  surf_normal_align_energy(const mati_t &surf, const matd_t &nods, const double eps, const double w)
      : surf_(surf), dim_(nods.size()), eps_(eps), w_(w) {
    // matd_t vert = zeros<double>(3, 3);
    // double sum_area = 0, value = 0;
    // for (size_t i = 0; i < surf.size(2); ++i) {
    //   vert = nods(colon(), surf(colon(), i));
    //   triangle_area_(&value, &vert[0]);
    //   sum_area += value;
    // }
    // w_ /= sum_area;
  }
  size_t Nx() const {
    return dim_;
  }
  int Val(const double *x, double *val) const {
    itr_matrix<const double *> X(3, dim_/3, x);
    matd_t vert = zeros<double>(3, 3); double value = 0;
    for (size_t i = 0; i < surf_.size(2); ++i) {
      vert = X(colon(), surf_(colon(), i));
      surf_normal_align_(&value, &vert[0], &eps_);
      *val += w_*value;
    }
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    itr_matrix<const double *> X(3, dim_/3, x);
    itr_matrix<double *> G(3, dim_/3, gra);
    matd_t vert = zeros<double>(3, 3), g = zeros<double>(3, 3);
    for (size_t i = 0; i < surf_.size(2); ++i) {
      vert = X(colon(), surf_(colon(), i));
      surf_normal_align_jac_(&g[0], &vert[0], &eps_);
      G(colon(), surf_(colon(), i)) += w_*g;
    }
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    itr_matrix<const double *> X(3, dim_/3, x);
    matd_t vert = zeros<double>(3, 3), H = zeros<double>(9, 9), Htmp = zeros<double>(9, 9);
    matd_t e = zeros<double>(9, 1), diag = zeros<double>(9, 9);
    for (size_t i = 0; i < surf_.size(2); ++i) {
      vert = X(colon(), surf_(colon(), i));
      for (int m = 0; m < 3; ++m) {
        area_normal_align_hes(&vert[0], eps_, m, &Htmp[0]);
        eig(Htmp, e);
        for (size_t k = 0; k < 9; ++k)
          diag(k, k) = std::max(0.0, e[k]);
        Htmp = temp(Htmp*temp(diag*trans(Htmp)));
        H += Htmp;
      }
      for (size_t p = 0; p < 9; ++p) {
        for (size_t q = 0; q < 9; ++q) {
          const size_t I = 3*surf_(p/3, i)+p%3;
          const size_t J = 3*surf_(q/3, i)+q%3;
          hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
        }
      }
    }
    return 0;
  }
  void scaling_args(const double scale_w, const double scale_eps) {
    w_ = w_*scale_w;
    eps_ = std::max(0.01, eps_*scale_eps);
  }
  double get_eps() const {
    return eps_;
  }
private:
  const mati_t &surf_;
  const size_t dim_;
  double eps_;
  double w_;
};

class tet_distortion_energy : public Functional<double>
{
public:
  tet_distortion_energy(const mati_t &tets, const matd_t &nods, const double w)
      : tets_(tets), w_(w), dim_(nods.size()) {
    vol_.resize(tets.size(2));
    Dm_.resize(9, tets.size(2));
    R_.resize(9, tets.size(2));
    #pragma omp parallel for
    for (size_t i = 0; i < tets.size(2); ++i) {
      matd_t base = nods(colon(), tets(colon(1, 3), i))-nods(colon(), tets(0, i))*ones<double>(1, 3);
      matd_t cp = base;
      vol_[i] = fabs(det(cp))/6.0;
      inv(base);
      Dm_(colon(), i) = base(colon());
    }
    double vol_sum = sum(vol_);
    vol_ /= vol_sum;
  }
  size_t Nx() const {
    return dim_;
  }
  int Val(const double *x, double *val) const {
    itr_matrix<const double *> X(3, dim_/3, x);

    matd_t vert = zeros<double>(3, 4); double value = 0;
    for (size_t i = 0; i < tets_.size(2); ++i) {
      vert = X(colon(), tets_(colon(), i));
      tet_distortion_(&value, &vert[0], &Dm_(0, i), &R_(0, i));
      *val += w_*vol_[i]*value;
    }
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    itr_matrix<const double *> X(3, dim_/3, x);
    itr_matrix<double *> G(3, dim_/3, gra);

    matd_t vert = zeros<double>(3, 4), g = zeros<double>(3, 4);
    for (size_t i = 0; i < tets_.size(2); ++i) {
      vert = X(colon(), tets_(colon(), i));
      tet_distortion_jac_(&g[0], &vert[0], &Dm_(0, i), &R_(0, i));
      G(colon(), tets_(colon(), i)) += w_*vol_[i]*g;
    }
    return 0;
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    matd_t H = zeros<double>(12, 12);
    for (size_t i = 0; i < tets_.size(2); ++i) {
      tet_distortion_hes_(&H[0], nullptr, &Dm_(0, i), nullptr);
      for (size_t p = 0; p < 12; ++p) {
        for (size_t q = 0; q < 12; ++q) {
          if ( H(p, q) != 0.0 ) {
            const size_t I = 3*tets_(p/3, i)+p%3;
            const size_t J = 3*tets_(q/3, i)+q%3;
            hes->push_back(Triplet<double>(I, J, w_*vol_[i]*H(p, q)));
          }
        }
      }
    }
    return 0;
  }
  void update_rotation(const double *x) {
    itr_matrix<const double *> X(3, dim_/3, x);
    #pragma omp parallel for
    for (size_t i = 0; i < tets_.size(2); ++i) {
      matd_t vert = X(colon(), tets_(colon(), i));
      matd_t F(3, 3), U(3, 3), S(3, 3), VT(3, 3), R(3, 3);
      calc_def_grad(&vert[0], &Dm_(0, i), &F[0]);
      svd(F, U, S, VT);
      R = U*VT;
      R_(colon(), i) = R(colon());
    }
  }
private:
  const mati_t &tets_;
  const size_t dim_;
  double w_;
  matd_t Dm_, vol_, R_;
};

//===============================================================================
//-------------------------------------------------------------------------------
//===============================================================================

int polycube_solver::unit_test() const {
  {
    matd_t nods = zeros<double>(3, 3);
    double area = 0;
    nods(0, 1) = 1; nods(1, 1) = -1;
    nods(0, 2) = 1; nods(1, 2) = 1;
    triangle_area_(&area, &nods[0]);
    cout << "area: " << area << endl;
  }
  {
    srand(time(NULL));
    matd_t nods = rand(3, 4);
    matd_t dm = nods(colon(), colon(1, 3))-nods(colon(), 0)*ones<double>(1, 3);
    inv(dm);

    matd_t F = (nods(colon(), colon(1, 3))-nods(colon(), 0)*ones<double>(1, 3))*dm;
    matd_t U, S, VT;
    svd(F, U, S, VT);
    matd_t R = U*VT;
    double value = 0;
    tet_distortion_(&value, &nods[0], &dm[0], &R[0]);
    cout << "distortion: " << value << endl;
  }
  {
    matd_t nods = rand(3, 3);
    double value = 0, eps = 0.1;
    surf_normal_align_(&value, &nods[0], &eps);
    cout << "normal align: " << value << endl;
  }
  {
    matd_t A = rand(3, 3);
    matd_t ATA = trans(A)*A;
    cout << ATA << endl;
    
    matd_t e(3, 1);
    eig(ATA, e);
    matd_t diag(3, 3);
    for (size_t i = 0; i < 3; ++i)
      diag(i, i) = e[i];
    cout << ATA*diag*trans(ATA) << endl;
  }
  {
    matd_t A = rand(3, 3);
    const double eps = 1;
    matd_t H(9, 9);
    surf_normal_align_hes_(&H[0], &A[0], &eps);
    cout << "H: " << H << endl;
    
    matd_t Ht = zeros<double>(9, 9);
    for (int i = 0; i < 3; ++i) {
      matd_t tmp = zeros<double>(9, 9);
      area_normal_align_hes(&A[0], eps, i, &tmp[0]);
      Ht += tmp;
    }
    cout << "Ht: " << Ht << endl;
  }
}

polycube_solver::polycube_solver(const mati_t &tets, const matd_t &nods, ptree &pt)
    : tets_(tets), pt_(pt) {
  shared_ptr<face2tet_adjacent> f2t(face2tet_adjacent::create(tets));
  bool check_order = true;
  jtf::mesh::get_outside_face(*f2t, surf_, check_order, &nods);

  const double eps = pt.get<double>("epsilon.value");
  const double w1 = pt.get<double>("weight.onenorm.value");
  const double wd = pt.get<double>("weight.distortion.value");
  
  buffer_.resize(2);
  buffer_[0] = alig_ = make_shared<surf_normal_align_energy>(tets, nods, eps, w1);
  buffer_[1] = arap_ = make_shared<tet_distortion_energy>(tets, nods, wd);
  energy_ = make_shared<energy_t<double>>(buffer_);
  
  area_cons_ = make_shared<area_preserving_cons>(surf_, nods);
}

double polycube_solver::eval_polycube_error(const matd_t &x) const {
  double area_sum = 0, align_sum = 0;
  double area = 0, align = 0;
  matd_t vert = zeros<double>(3, 3);
  const double eps = alig_->get_eps();
  
  for (size_t i = 0; i < surf_.size(2); ++i) {
    vert = x(colon(), surf_(colon(), i));
    surf_normal_align_(&align, &vert[0], &eps);
    triangle_area_(&area, &vert[0]);
    area_sum += area;
    align_sum += align;
  }

  return align_sum-area_sum;
}

int polycube_solver::deform(matd_t &x) const {
  const size_t dim = energy_->Nx();
  ASSERT(x.size() == dim);
  Map<VectorXd> xstar(&x[0], dim);

  CholmodSimplicialLDLT<SparseMatrix<double>> solver;

  const size_t maxits = pt_.get<size_t>("maxits.value");
  matd_t error = zeros<double>(maxits, 1);
  size_t args_update_count = 0, d0 = 0, count = 0;
  for (size_t iter = 0; iter < maxits; ++iter)
  {
    // evaluate polycube error
    error[iter] = eval_polycube_error(x);
    if ( iter % 10 == 0 )
      cout << "\t # iter " << iter << ", polycube error: " << error[iter] << endl;

    // convergence test
    if ( fabs(error[iter]) < 0.001 || args_update_count > 20 ) {
      cout << "[INFO] polycube solver converged after " << iter << " iterations" << endl;
      break;
    }
    
    // update parameters for L1 energy
    if ( count == 20 ) //>= 5 && sum(fabs(error(colon(iter-4, iter))-error(colon(iter-5, iter-1)))) < fabs(error[d0])/100 ) {
    {
      cout << "\t ! UPDATE" << endl;
      ++args_update_count;
      count = 0;
      d0 = iter+1;
      alig_->scaling_args(2, 0.5);
    }
    ++count;
    
    // update rotation
    arap_->update_rotation(&x[0]);

    // Schur complement for constrained optimization
    double ve = 0; {
      energy_->Val(&x[0], &ve);
      if ( iter % 10 == 0 ) {
        cout << "\t\t energy value: " << ve << endl;
        double vtmp = 0;
        alig_->Val(&x[0], &vtmp);
        cout << "\t\t alignment: " << vtmp << endl;
        cout << "\t\t distortion: " << ve-vtmp << endl;
      }
    }
    VectorXd g = VectorXd::Zero(dim); {
      energy_->Gra(&x[0], g.data());
    }
    SparseMatrix<double> H(dim, dim); {
      vector<Triplet<double>> trips;
      energy_->Hes(&x[0], &trips);
      H.reserve(trips.size());
      H.setFromTriplets(trips.begin(), trips.end());
    }
    double vc = 0; {
      area_cons_->Val(&x[0], &vc);
      if ( iter % 10 == 0 )
        cout << "\t\t constraint value: " << vc << endl << endl;
    }
    VectorXd gc = VectorXd::Zero(dim); {
      area_cons_->Gra(&x[0], gc.data());
    }

    solver.compute(H);
    ASSERT(solver.info() == Success);

    double lambda = (-gc.dot(solver.solve(g))+vc)/(gc.dot(solver.solve(gc)));
    VectorXd dx = -solver.solve(g)-lambda*solver.solve(gc);
    xstar += dx;
  }
  
  return 0;
}

}
