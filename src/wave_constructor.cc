#include "wave_constructor.h"

#include <iostream>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>

#include "def.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;
using namespace jtf::mesh;
using mati_t = zjucad::matrix::matrix<size_t>;
using matd_t = zjucad::matrix::matrix<double>;

namespace riemann {

extern "C" {

void wave_value_condition_(double *val, const double *f, const double *c);
void wave_value_condition_jac_(double *jac, const double *f, const double *c);
void wave_value_condition_hes_(double *hes, const double *f, const double *c);

void modulus_condition_(double *val, const double *f);
void modulus_condition_jac_(double *jac, const double *f);
void modulus_condition_hes_(double *hes, const double *f);

void phase_condition_(double *val, const double *f);
void phase_condition_jac_(double *jac, const double *f);
void phase_condition_hes_(double *hes, const double *f);

}

class wave_value_energy : public Functional<double>
{
public:
  wave_value_energy(const mati_t &edge, const matd_t &f, const matd_t &c, const double w=1.0)
    : edge_(edge), dim_(f.size()), c_(c), w_(w) {}
  size_t Nx() const {
    return dim_;
  }
  int Val(const double *x, double *val) const {
    itr_matrix<const double *> X(4, dim_/4, x);
    for (size_t i = 0; i < edge_.size(2); ++i) {
      matd_t vert = X(colon(), edge_(colon(), i));
      double value = 0;
      wave_value_condition_(&value, &vert[0], &c_(0, i));
      *val += w_*value;
    }
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    itr_matrix<const double *> X(4, dim_/4, x);
    itr_matrix<double *> G(4, dim_/4, gra);
    for (size_t i = 0; i <edge_.size(2); ++i) {
      matd_t vert = X(colon(), edge_(colon(), i));
      matd_t g = zeros<double>(4, 2);
      wave_value_condition_jac_(&g[0], &vert[0], &c_(0, i));
      G(colon(), edge_(colon(), i)) += w_*g;
    }
    return 0;
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    itr_matrix<const double *> X(4, dim_/4, x);
    for (size_t i = 0; i < edge_.size(2); ++i) {
      matd_t vert = X(colon(), edge_(colon(), i));
      matd_t H = zeros<double>(8, 8);
      wave_value_condition_hes_(&H[0], &vert[0], &c_(0, i));
      for (size_t q = 0; q < H.size(2); ++q) {
        for (size_t p = 0; p < H.size(1); ++p) {
          if ( H(p, q) != 0.0 ) {
            size_t I = 4*edge_(p/4, i)+p%4;
            size_t J = 4*edge_(q/4, i)+q%4;
            hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
          }
        }
      }
    }
    return 0;
  }
private:
  const mati_t &edge_;
  const matd_t &c_;
  const size_t dim_;
  const double w_;
};

class modulus_energy : public Functional<double>
{
public:
  modulus_energy(const matd_t &f, const double w=0.15)
    : dim_(f.size()), w_(w) {}
  size_t Nx() const {
    return dim_;
  }
  int Val(const double *x, double *val) const {
    itr_matrix<const double *> X(4, dim_/4, x);
    for (size_t i = 0; i < X.size(2); ++i) {
      double value = 0;
      modulus_condition_(&value, &X(0, i));
      *val += w_*value;
    }
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    itr_matrix<const double *> X(4, dim_/4, x);
    itr_matrix<double *> G(4, dim_/4, gra);
#pragma omp parallel for
    for (size_t i = 0; i < X.size(2); ++i) {
      matd_t g  = zeros<double>(4, 1);
      modulus_condition_jac_(&g[0], &X(0, i));
      G(colon(), i) += w_*g;
    }
    return 0;
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    itr_matrix<const double *> X(4, dim_/4, x);
    for (size_t i = 0; i < X.size(2); ++i) {
      matd_t H = zeros<double>(4, 4);
      modulus_condition_hes_(&H[0], &X(0, i));
      for (size_t q = 0; q < H.size(2); ++q) {
        for (size_t p = 0; p < H.size(1); ++p)  {
          if ( H(p, q) != 0.0 )
            hes->push_back(Triplet<double>(4*i+p, 4*i+q, w_*H(p, q)));
        }
      }
    }
    return 0;
  }
private:
  const size_t dim_;
  const double w_;
};

class phase_energy : public Functional<double>
{
public:
  phase_energy(const matd_t &f, const double w=0.15)
    : dim_(f.size()), w_(w) {}
  size_t Nx() const {
    return dim_;
  }
  int Val(const double *x, double *val) const {
    itr_matrix<const double *> X(4, dim_/4, x);
    for (size_t i = 0; i < X.size(2); ++i) {
      double value = 0;
      phase_condition_(&value, &X(0, i));
      *val += w_*value;
    }
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    itr_matrix<const double *> X(4, dim_/4, x);
    itr_matrix<double *> G(4, dim_/4, gra);
#pragma omp parallel for
    for (size_t i = 0; i < X.size(2); ++i) {
      matd_t g  = zeros<double>(4, 1);
      phase_condition_jac_(&g[0], &X(0, i));
      G(colon(), i) += w_*g;
    }
    return 0;
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    itr_matrix<const double *> X(4, dim_/4, x);
    for (size_t i = 0; i < X.size(2); ++i) {
      matd_t H = zeros<double>(4, 4);
      phase_condition_hes_(&H[0], &X(0, i));
      for (size_t q = 0; q < H.size(2); ++q) {
        for (size_t p = 0; p < H.size(1); ++p)  {
          if ( H(p, q) != 0.0 )
            hes->push_back(Triplet<double>(4*i+p, 4*i+q, w_*H(p, q)));
        }
      }
    }
    return 0;
  }
private:
  const size_t dim_;
  const double w_;
};
//==============================================================================
/// @todo CHECK!
int wave_constructor::solve_phase_transition() {
  alpha_ij_.resize(1,  edge_.size(2));
  beta_ij_.resize(1, edge_.size(2));
  for (size_t i = 0; i < edge_.size(2); ++i) {
    const size_t I = edge_(0, i);
    const size_t J = edge_(1, i);
    alpha_ij_[i] = w*dot(nods_(colon(), I)-nods_(colon(), J), frameX);
    beta_ij_[i] = w*dot(nods_(colon(), I)-nods_(colon(), J), frameY);
  }
  c_ij_.resize(4, edge_.size(2));
  for (size_t i = 0; i < c_ij_.size(2); ++i) {
    c_ij_(0, i) = cos(alpha_ij_[i])*cos(beta_ij_[i]);
    c_ij_(1, i) = -cos(alpha_ij_[i])*sin(beta_ij_[i]);
    c_ij_(2, i) = -sin(alpha_ij_[i])*cos(beta_ij_[i]);
    c_ij_(3, i) = sin(alpha_ij_[i])*sin(beta_ij_[i]);
  }
  return 0;
}

//==============================================================================
void wave_constructor::test_wave_conditions() const {
  const double f[4] = {1, 2, 3, 4};
  double val = 0;
  modulus_condition_(&val, f);
  cout << val << endl;
  cout << pow(2*2+3*3+4*4, 2) << endl;
}

}
