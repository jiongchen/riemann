#include "advanced_mips.h"

#include <iostream>
#include <unordered_map>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <zjucad/matrix/itr_matrix.h>
#include <jtflib/mesh/mesh.h>
#include <boost/math/special_functions/sign.hpp>

#include "config.h"
#include "geometry_extend.h"

using namespace std;
using namespace zjucad::matrix;

namespace riemann {

extern "C" {

void mips_2d_(double *val, const double *x, const double *D);
void mips_2d_jac_(double *jac, const double *x, const double *D);
void mips_2d_hes_(double *hes, const double *x, const double *D);

void det_2d_(double *val, const double *x, const double *D);
void det_2d_jac_(double *jac, const double *x, const double *D);
void det_2d_hes_(double *hes, const double *x, const double *D);

void advanced_iso_2d_(double *val, const double *x, const double *D, const double *s);
void advanced_iso_2d_jac_(double *jac, const double *x, const double *D, const double *s);
void advanced_iso_2d_hes_(double *hes, const double *x, const double *D, const double *s);

}

class mips_energy
{
public:
  virtual size_t dim() const = 0;
  virtual int val(const double *x, double *value) const = 0;
  virtual int val(const double *x, const size_t id, double *f) const = 0;
  virtual int gra(const double *x, const size_t id, double *g) const = 0;
  virtual int hes(const double *x, const size_t id, double *H) const = 0;
};

class mips_energy_2d : public mips_energy
{
public:
  mips_energy_2d(const mati_t &tris, const matd_t &nods)
    : dim_(nods.size()), tris_(tris), s_(5.0) {
    binv_.resize(4, tris_.size(2));
    for (size_t i = 0; i < tris.size(2); ++i) {
      matd_t base = nods(colon(), tris_(colon(1, 2), i))-nods(colon(), tris_(0, i))*ones<double>(1, 2);
      inv(base);
      std::copy(base.begin(), base.end(), &binv_(0, i));
    }
    calc_one_ring_face(tris_, p2f_);
  }
  size_t dim() const {
    return dim_;
  }
  int val(const double *x, double *value) const {
    itr_matrix<const double *> X(2, dim_/2, x);
    for (size_t i = 0; i < tris_.size(2); ++i) {
      matd_t vert = X(colon(), tris_(colon(), i));
      double va = 0;
      advanced_iso_2d_(&va, &vert[0], &binv_(0, i), &s_);
      *value += va;
    }
    return 0;
  }
  int val(const double *x, const size_t id, double *v) const {
    itr_matrix<const double *> X(2, dim_/2, x);
    for (auto &fa : p2f_[id]) {
      matd_t vert = X(colon(), tris_(colon(), fa));
      double va = 0;
      advanced_iso_2d_(&va, &vert[0], &binv_(0, fa), &s_);
      *v += va;
    }
    return 0;
  }
  int gra(const double *x, const size_t id, double *g) const {
    itr_matrix<const double *> X(2, dim_/2, x);
    itr_matrix<double *> grad(2, 1, g);
    for (auto &fa : p2f_[id]) {
      matd_t vert = X(colon(), tris_(colon(), fa));
      matd_t jac = zeros<double>(2, 3);
      advanced_iso_2d_jac_(&jac[0], &vert[0], &binv_(0, fa), &s_);
      int pos = std::find(&tris_(0, fa), &tris_(0, fa)+3, id)-&tris_(0, fa);
      ASSERT(pos >= 0 && pos < 3);
      grad += jac(colon(), pos);
    }
    return 0;
  }
  int hes(const double *x, const size_t id, double *H) const {
    return __LINE__;
  }
public:
  vector<vector<size_t>> p2f_;
private:
  const size_t dim_;
  const mati_t &tris_;
  const double s_;
  matd_t binv_;
};

class move_vertex
{
public:
  move_vertex(const mati_t &tris, const matd_t &nods)
    : tris_(tris) {
    energy_ = make_shared<mips_energy_2d>(tris, nods);
    matd_t e0 = nods(colon(), tris(colon(1, 2), 0))-nods(colon(), tris(colon(0, 1), 0));
    sign_ = boost::math::sign(det(e0));
  }
  int operator()(const size_t id, double *x) const {
    itr_matrix<double *> X(2, energy_->dim()/2, x);
    double v0 = 0, v1 = 0;
    energy_->val(x, id, &v0);
    matd_t g = zeros<double>(2, 1);
    energy_->gra(x, id, &g[0]);
    double lambda = calc_step_length(id, x, g);
//    cout << "\t\t#pid: " << id << " step: " << lambda << endl;
    X(colon(), id) -= lambda*g;
    energy_->val(x, id, &v1);
    if ( v1 > v0 || has_invert_elem_one_ring(id, x) )
      X(colon(), id) += 0.85*lambda*g;
    return 0;
  }
public:
  shared_ptr<mips_energy_2d> energy_;
private:
  bool has_invert_elem_one_ring(const size_t id, const double *x) const {
    itr_matrix<const double *> X(2, energy_->dim()/2, x);
    for (auto &fa : energy_->p2f_[id]) {
      matd_t ev = X(colon(), tris_(colon(1, 2), fa))-X(colon(), tris_(colon(0, 1), fa));
      if ( det(ev)*sign_ <= 0.0 )
        return true;
    }
    return false;
  }
  double calc_step_length(const size_t id, const double *x, matd_t &d) const {
    itr_matrix<const double *> X(2, energy_->dim()/2, x);
    for (auto &fa : energy_->p2f_[id]) {
      int pos = std::find(&tris_(0, fa), &tris_(0, fa)+3, id)-&tris_(0, fa);
      ASSERT(pos >= 0 && pos < 3);
      matd_t e0 = X(colon(), tris_((pos+1)%3, fa))-X(colon(), tris_(pos, fa));
      matd_t e1 = X(colon(), tris_((pos+2)%3, fa))-X(colon(), tris_(pos, fa));
      if ( dot(cross(e0, d), cross(d, e1)) > 0 ) {
        return (e0(0, 0)*e1(1, 0)-e1(0, 0)*e0(1, 0))
            /(d(0, 0)*(e1(1, 0)-e0(1, 0))+d(1, 0)*(e0(0, 0)-e1(0, 0)));
      }
    }
    return 1.0;
  }
  const mati_t &tris_;
  double sign_;
};

//========================== 2d mesh deformer ==================================
mips_deformer_2d::mips_deformer_2d(const mati_t &tris, const matd_t &nods)
  : dim_(nods.size()) {
  move_ = make_shared<move_vertex>(tris, nods);
}

void mips_deformer_2d::set_fixed_vert(const unordered_set<size_t> &fixed) {
  fixed_ = fixed;
}

int mips_deformer_2d::deform(double *x, const size_t maxiter) const {
  cout << "[info] deform the mesh" << endl;
  double prev_val = 0, post_val = 0;
  move_->energy_->val(x, &prev_val);
  // iterative solve
  for (size_t iter = 0; iter < maxiter; ++iter) {
    if ( iter % 1 == 0 ) {
      cout << "\t@iter " << iter << " energy: " << prev_val << endl;
    }
    apply(x);
    move_->energy_->val(x, &post_val);
    if ( std::fabs(post_val-prev_val) < 1e-6*std::fabs(prev_val) ) {
      cout << "\t@converged after " << iter+1 << " iteration\n";
      break;
    }
    prev_val = post_val;
  }
  return 0;
}

int mips_deformer_2d::apply(double *x) const {
  for (size_t i = 0; i < dim_/2; ++i) {
    if ( fixed_.find(i) == fixed_.end() )
      (*move_)(i, x);
  }
  return 0;
}

//==============================================================================
void mips_deformer_2d::unit_test() const {
  matd_t nods = rand(2, 3);
  matd_t base = nods(colon(), colon(1, 2))-nods(colon(), 0)*ones<double>(1, 2);
  inv(base);

  double angle = M_PI/3;
  matd_t R(2, 2);
  R(0, 0) = cos(angle);
  R(0, 1) = -sin(angle);
  R(1, 0) = sin(angle);
  R(1, 1) = cos(angle);
  nods = temp(R*nods);
  nods *= 1;

  double value = 0;
  mips_2d_(&value, &nods[0], &base[0]);
  cout << value << endl;

  det_2d_(&value, &nods[0], &base[0]);
  cout << value << endl;

  double s = 5;
  advanced_iso_2d_(&value, &nods[0], &base[0], &s);
  cout << value << endl;
}

}
