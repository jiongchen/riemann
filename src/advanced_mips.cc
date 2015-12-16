#include "advanced_mips.h"

#include <iostream>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <zjucad/matrix/itr_matrix.h>

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

}

class mips_energy
{
public:
  virtual int val(const double *x, double *value) const = 0;
  virtual int gra(const double *x, const size_t id, double *grad) const = 0;
  virtual int hes(const double *x, const size_t id, double *hess) const = 0;
};

static void gen_rel_elem_to_vert(const mati_t &cell, vector<vector<size_t>> &result) {
  const size_t vert_num = max(cell)+1;
  result.resize(vert_num);
  for (size_t j = 0; j < cell.size(2); ++j) {
    for (size_t i = 0; i < cell.size(1); ++i) {
      result[cell(i, j)].push_back(j);
    }
  }
}

class mips_energy_2d : public mips_energy
{
public:
  mips_energy_2d(const mati_t &tris, const matd_t &nods, const double w=1.0)
    : w_(w), dim_(nods.size()), tris_(tris) {
    binv_.resize(4, tris.size(2));
#pragma omp parallel for
    for (size_t i = 0; i < tris.size(2); ++i) {
      matd_t base = nods(colon(), colon(1, 2))-nods(colon(), 0)*ones<double>(1, 2);
      inv(base);
      copy(base.begin(), base.end(), &binv_(0, i));
    }
//    build_related_elem_map(tris, rel_elem_);
  }
  size_t dim() const {
    return dim_;
  }
  int val(const double *x, double *value) const {
    itr_matrix<const double *> X(2, dim_/2, x);
    for (size_t i = 0; i < tris_.size(2); ++i) {
      matd_t vert = X(colon(), tris_(colon(), i));
      double v0 = 0, v1 = 0;
      mips_2d_(&v0, &vert[0], &binv_(0, i));
      det_2d_(&v1, &vert[0], &binv_(0, i));
      *value += w_*(v0+v1);
    }
    return 0;
  }
  int gra(const double *x, const size_t id, double *grad) const {
    itr_matrix<const double *> X(2, dim_/2, x);
    itr_matrix<double *> G(2, 1, grad);
//    for (auto i : related_elem[id]) {
//      matd_t vert = X(colon(), tris_(colon(), i));
//      matd_t g = zeros<double>(6, 1);
//      mips_2d_jac_(&g[0], &vert[0], &binv(0, i));
//      std::find(tris(0, i), tris(3, i), id)-tris(0, );
//      grad += ;
//    }
    return 0;
  }
  int hes(const double *x, const size_t id, double *hess) const {
    itr_matrix<const double *> X(3, dim_/3, x);
    itr_matrix<double *> H(2, 2, hess);

    return 0;
  }
private:
  const double w_;
  const size_t dim_;
  const mati_t &tris_;
  matd_t binv_;
  vector<vector<size_t>> rel_elem_;
};

class move_vertex
{
public:
  move_vertex() {
//    energy = make_shared<>();
  }
  int operator()(const size_t id, double *x) {

  }
private:
  shared_ptr<mips_energy_2d> energy;
};
//========================== 2d mesh deformer ==================================
mips_deformer_2d::mips_deformer_2d(const mati_t &tris, const matd_t &nods) {

}



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
  nods *= 5;

  double value = 0;
  mips_2d_(&value, &nods[0], &base[0]);
  cout << value << endl;

  det_2d_(&value, &nods[0], &base[0]);
  cout << value << endl;
}

}
