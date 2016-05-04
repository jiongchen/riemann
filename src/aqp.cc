#include "aqp.h"

#include "def.h"

using namespace std;
using namespace Eigen;

namespace riemann {

extern "C" {

void two_dim_arap_(double *val, const double *x, const double *Dm, const double *R, const double *area);
void two_dim_arap_jac_(double *jac, const double *x, const double *Dm, const double *R, const double *area);
void two_dim_arap_hes_(double *hes, const double *x, const double *Dm, const double *R, const double *area);

void two_dim_iso_(double *val, const double *x, const double *Dm, const double *area);
void two_dim_iso_jac_(double *jac, const double *x, const double *Dm, const double *area);
void two_dim_iso_hes_(double *hes, const double *x, const double *Dm, const double *area);

}

class two_dim_arap_energy : public Functional<double>
{
public:
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, vector<Triplet<double>> *hes) const;
};

class two_dim_iso_energy : public Functional<double>
{
public:
  size_t Nx() const;
  int Val(const double *x, double *val) const;
  int Gra(const double *x, double *gra) const;
  int Hes(const double *x, vector<Triplet<double>> *hes) const;
};

class two_dim_pos_cons : public Constraint<double>
{
public:
  size_t Nx() const;
  size_t Nf() const;
  int Val(const double *x, double *val) const;
  int Jac(const double *x, const size_t off, vector<Triplet<double>> *jac) const;
};

}
