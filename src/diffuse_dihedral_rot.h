#ifndef DIFFUSE_DIHEDRAL_ROT_H
#define DIFFUSE_DIHEDRAL_ROT_H

#include <unordered_set>
#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>

namespace riemann {

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

class diffuse_arap_energy;
class graph_t;
typedef graph_t tree_t;

class diffuse_arap_encoder
{
public:
  void calc_delta_angle(const mati_t &tris, const matd_t &prev, const matd_t &curr,
                        const tree_t &g, const size_t root_face, matd_t &root_curr, std::vector<double> &da);
};

class diffuse_arap_decoder
{
public:
  diffuse_arap_decoder(const mati_t &tris, const matd_t &nods);
  int estimate_rotation(const matd_t &prev, const tree_t &g, const size_t root_face, const matd_t &root_nods, const std::vector<double> &da);
  int pin_down_vert(const size_t id, const double *pos);
  int solve(matd_t &curr);
private:
  const mati_t &tris_;
  const matd_t &nods_;
  size_t dim_;
  std::unordered_set<size_t> fixDoF_;
  std::shared_ptr<diffuse_arap_energy> energy_;
  std::vector<Eigen::Matrix3d> R_;
  Eigen::VectorXd X_;
};

}

#endif
