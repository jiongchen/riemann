#include <iostream>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "src/geometry_extend.h"
#include "src/write_vtk.h"
#include "src/lscm.h"
#include "src/arap_param.h"

using namespace std;
using namespace zjucad::matrix;
using namespace riemann;
using jtf::mesh::edge2cell_adjacent;

int main(int argc, char *argv[])
{
  boost::filesystem::create_directories("./arap_param");

  mati_t tris; matd_t nods, uv, uv3;
  jtf::mesh::load_obj(argv[1], tris, nods);
  uv.resize(2, nods.size(2));
  uv3.resize(3, nods.size(2));

  mati_t bnd;
  shared_ptr<edge2cell_adjacent> e2c(edge2cell_adjacent::create(tris, false));
  jtf::mesh::get_boundary_edge(*e2c, bnd);
  std::sort(bnd.begin(), bnd.end());

  // LSCM
  lscm_param lscm(tris, nods);
  const double pos[4] = {0, 0, 1, 0};
  lscm.set_fixed_bnd_vert(bnd[0], &pos[0]);
  lscm.set_fixed_bnd_vert(bnd[bnd.size()/2], &pos[2]);
  cout << "[INFO] fixed point: " << bnd[0] << " " << bnd[bnd.size()/2] << endl;
  lscm.apply(&uv[0]);
  uv3(colon(0, 1), colon()) = uv;
  jtf::mesh::save_obj("./arap_param/lscm_param.obj", tris, uv3);

  // ARAP
  arap_param_solver arap(tris, nods);
  arap.precompute();
  arap.solve(&uv[0]);
  uv3(colon(0, 1), colon()) = uv;
  jtf::mesh::save_obj("./arap_param/arap_param.obj", tris, uv3);

  cout << "[INFO] done\n";
  return 0;
}
