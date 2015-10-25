#include "shell.h"

#include <iostream>
#include <zjucad/matrix/itr_matrix.h>

#include "def.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

namespace riemann {

extern "C" {

void calc_edge_length_(double *val, const double *x);
void calc_edge_length_jac_(double *jac, const double *x);

void calc_dihedral_angle_(double *val, const double *x);
void calc_dihedral_angle_jac_(double *jac, const double *x);

void calc_volume_(double *val, const double *x);
void calc_volume_jac_(double *jac, const double *x);

}

class stretch_constraint : public Constraint<double>
{
public:
  stretch_constraint(const mati_t &edge, const matd_t &nods, const double w)
    : dim_(nods.size()), edges_(edge), w_(sqrt(w)) {
    len_.resize(edges_.size(2), 1);
  }
  size_t Nx() const {
    return dim_;
  }
  size_t Nf() const {
    return edges_.size(2);
  }
  int Val(const double *x, double *val) const {
    itr_matrix<const double *> X(3, Nx()/3, x);
    itr_matrix<double *> V(Nf(), 1, val);
    for (size_t i = 0; i < edges_.size(2); ++i) {

    }
    return 0;
  }
  int Jac(const double *x, const size_t off, vector<Triplet<double>> *jac) const {

  }
private:
  const size_t dim_;
  const double w_;
  const mati_t &edges_;
  matd_t len_;
};

class bending_constraint : public Constraint<double>
{
public:
  bending_constraint(const mati_t &diamond, const matd_t &nods, const double w);
  size_t Nx() const;
  size_t Nf() const;
  int Val(const double *x, double *val) const;
  int Jac(const double *x, const size_t off, vector<Triplet<double>> *jac) const;
};

class volume_constraint : public Constraint<double>
{
public:
};
//==============================================================================
void shell_solver::temp_test() const {
  const double x[6] = {0,0,0, 3,4,5};
  double len = 0;
  calc_edge_length_(&len, x);
  cout << "len^2: " << len*len << endl;

  const double y[12] = {1,0,0, 0,0,0, 0,1,0, 0,0,-11};
  double angle = 0;
  calc_dihedral_angle_(&angle, y);
  cout << "angle: " << 180-angle/M_PI*180 << endl;
}

}
