#ifndef PARAM_LSCM_H
#define PARAM_LSCM_H

#include <zjucad/matrix/matrix.h>
#include <unordered_set>

#include "def.h"

namespace riemann {

class lscm_param
{
public:
  typedef zjucad::matrix::matrix<size_t> mati_t;
  typedef zjucad::matrix::matrix<double> matd_t;
  lscm_param(const mati_t &tris, const matd_t &nods);
  ///> used for remove the kernel of Lc
  ///> in spetral conformal map, dont need to pinn vertices on boundary
  void set_fixed_bnd_vert(const size_t id, const double *x);
  int apply();
  int apply_spetral();
  int get_param_mesh(mati_t *param_tris, matd_t *param_nods);
private:
  int get_boundary_loop(const mati_t &cell, std::vector<size_t> &bnd);
  std::vector<std::shared_ptr<Functional<double>>> buff_;
  std::shared_ptr<Functional<double>> conformal_energy_;
  const mati_t &tris_;
  std::vector<size_t> bnd_;
  std::unordered_set<size_t> fixed_dofs_;
  std::vector<size_t> g2l_;
  Eigen::VectorXd uv_;
};

}
#endif
