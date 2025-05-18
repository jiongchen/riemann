#include <iostream>
#include <boost/filesystem.hpp>
#include <jtflib/mesh/io.h>
#include <zjucad/matrix/matrix.h>

#include "../src/deform_transfer.h"

using namespace std;
using namespace riemann;

int main(int argc, char *argv[])
{
  if ( argc != 5 ) {
    cerr << "usage: " << argv[0] << " source_ref.obj target_ref.obj vert_makers.cons source_def.obj\n";
    return __LINE__;
  }
  boost::filesystem::create_directories("./result/dt");

  deform_transfer dt;
  dt.load_reference_source_mesh(argv[1]);
  dt.load_reference_target_mesh(argv[2]);
  dt.load_vertex_markers(argv[3]);

  dt.init();

  dt.see_source_markers("./result/dt/src_markers.vtk");
  dt.see_target_markers("./result/dt/tar_markers.vtk");
  dt.see_ghost_tet_mesh("./result/dt/ghost_sr.vtk", "source_ref");
  dt.see_ghost_tet_mesh("./result/dt/ghost_tr.vtk", "target_ref");
  dt.save_reference_source_mesh("./result/dt/source_ref.obj");
  dt.save_reference_target_mesh("./result/dt/target_ref.obj");

//  dt.solve_corres_harmonic();

  dt.solve_corres_precompute();
  dt.solve_corres_first_phase();
  dt.see_corres_mesh("./result/dt/out_first.obj");
  dt.solve_corres_second_phase();
  dt.see_corres_mesh("./result/dt/out_second.obj");
  dt.compute_triangle_corres();

  dt.deformation_transfer_precompute();

  dt.load_deformed_source_mesh(argv[4]);
  dt.deformation_transfer();
  dt.save_deformed_source_mesh("./result/dt/source_def.obj");
  dt.save_deformed_target_mesh("./result/dt/target_def.obj");

  cout << "[info] all done\n";
  return 0;
}
