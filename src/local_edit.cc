#include "local_edit.h"

#include <jtflib/mesh/mesh.h>

#include "def.h"
#include "util.h"

using namespace std;
using namespace Eigen;
using namespace zjucad::matrix;

namespace riemann {

constrained_mesh_editor::constrained_mesh_editor(const mati_t &quad, const matd_t &nods)
  : quad_(quad), nods_(nods) {

}

int constrained_mesh_editor::set_handles(const Json::Value &json) {
  if ( json.empty() ) {
    cerr << "[INFO] no handles\n";
    return __LINE__;
  }
  ASSERT(json.isArray());
  for (int i = 0; i < json.size(); ++i) {
    size_t vid = json[i]["vid"].asInt();
    ASSERT(json[i]["disp"].isArray());
    double dx = json[i]["disp"][0].asDouble();
    double dy = json[i]["disp"][1].asDouble();
    double dz = json[i]["disp"][2].asDouble();
    handle_.insert(make_pair(vid, Vector3d(dx, dy, dz)));
  }
  return 0;
}

int constrained_mesh_editor::deform(double *d, const Json::Value &json) const {
  return 0;
}

}
