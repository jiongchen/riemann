#include <iostream>
#include <boost/filesystem.hpp>

#include "src/frame_field_deform.h"

using namespace std;
using namespace geom_deform;

int main(int argc, char *argv[])
{
  if ( argc != 2 ) {
    cerr << "# usage: " << argv[0] << " model.obj\n";
    return __LINE__;
  }
  boost::filesystem::create_directory("./ff_deform");

  frame_field_deform deformer;
  deformer.load_mesh(argv[1]);
  deformer.save_original_mesh("./ff_deform/origin.obj");
  deformer.save_local_frame("./ff_deform/local_frame.vtk", 0.005);

  cout << "[INFO] done\n";
  return 0;
}
