#ifndef DIFFUSE_DIHEDRAL_ROT_H
#define DIFFUSE_DIHEDRAL_ROT_H

#include <unordered_set>
#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>

namespace riemann {

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

class graph_t;
typedef graph_t tree_t;

void diffuse_rotation(const mati_t &tris, const matd_t &vrest, const matd_t &vcurr,
                      const size_t root_face, const tree_t &g, std::vector<Eigen::Matrix3d> &rot);

class diffuse_arap_energy;

class diffuse_arap_solver
{
public:
  diffuse_arap_solver(const mati_t &tets, const matd_t &nods, const std::vector<Eigen::Matrix3d> &R);
  int pin_down_vert(const size_t id, const double *pos, double *x);
  int solve(double *x) const;
private:
  const size_t dim_;
  std::unordered_set<size_t> fixDoF_;
  std::shared_ptr<diffuse_arap_energy> energy_;
};

}

#endif
