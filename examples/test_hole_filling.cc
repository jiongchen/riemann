#include <iostream>
#include <boost/filesystem.hpp>

#include "igl/readOBJ.h"
#include "igl/writeOBJ.h"
#include "igl/boundary_loop.h"

#include "src/vtk.h"
#include "src/hole_filling.h"
#include "src/vtk.h"
#include "src/write_vtk.h"

using namespace std;
using namespace Eigen;
using namespace riemann;

int main(int argc, char *argv[])
{
  if ( argc != 2 ) {
    cerr << "# usage: prog model.obj\n";
    return __LINE__;
  }
  boost::filesystem::create_directory("./hole_filling");

  mati tris; matd nods;
  igl::readOBJ(argv[1], nods, tris);
  igl::writeOBJ("./hole_filling/origin.obj", nods, tris);

  boundary_loop loop;
  get_boundary_loop(tris, nods, loop);
  calc_bnd_local_frame(tris, nods, loop);

  loop.T *= 0.1;
  loop.B *= 0.1;
  loop.N *= 0.1;
  draw_vert_direct_field("./hole_filling/T.vtk", loop.pos.data(), loop.pos.rows(), loop.T.data());
  draw_vert_direct_field("./hole_filling/B.vtk", loop.pos.data(), loop.pos.rows(), loop.B.data());
  draw_vert_direct_field("./hole_filling/N.vtk", loop.pos.data(), loop.pos.rows(), loop.N.data());

  cout << "[INFO] done\n";
  return 0;
}
