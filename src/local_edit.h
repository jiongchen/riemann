#ifndef LOCAL_MODIFICATIONS
#define LOCAL_MODIFICATIONS

#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>
#include <unordered_map>

#include "json.h"

namespace riemann {

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

/**
 * @brief reference: Exploring Local Modifications for Constrained Meshes
 */
class constrained_mesh_editor
{
public:
  constrained_mesh_editor(const mati_t &tris, const matd_t &nods);
  int set_handles(const Json::Value &json);
  int deform(double *d) const;
private:
  std::unordered_map<size_t, Eigen::Vector3d> handle_;
};

}

#endif
