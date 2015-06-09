#include <iostream>
#include <boost/filesystem.hpp>

#include "src/frame_field_deform.h"

using namespace std;
using namespace geom_deform;

int main(int argc, char *argv[])
{
  if ( argc != 3 ) {
    cerr << "# usage: " << argv[0] << " model.obj cons.ff\n";
    return __LINE__;
  }
  boost::filesystem::create_directory("./ff_deform");

  frame_field_deform deformer;
  deformer.load_mesh(argv[1]);
  deformer.save_original_mesh("./ff_deform/origin.obj");
  deformer.visualize_local_bases("./ff_deform/local_frame.vtk", 0.005);

  deformer.load_constraints(argv[2]);
  deformer.visualize_init_frames("./ff_deform/init_ff.vtk");
  deformer.visualize_frame_fields("./ff_deform/interp_ff.vtk");

  cout << "[INFO] done\n";
  return 0;
}
