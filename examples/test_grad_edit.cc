#include <iostream>
#include <boost/filesystem.hpp>

#include "src/gradient_deform.h"

using namespace std;
using namespace geom_deform;

int main(int argc, char *argv[])
{
  if ( argc != 2 ) {
    cerr << "# Usage: " << argv[0] << " model.obj\n";
    return __LINE__;
  }
  boost::filesystem::create_directory("./grad_edit");

  gradient_field_deform handle;
  handle.load_origin_model(argv[1]);
  handle.save_origin_model("./grad_edit/origin.obj");
  handle.init();
  enum {X, Y, Z};
  handle.see_coord_grad_fields("./grad_edit/grad_x.vtk", X);
  handle.see_coord_grad_fields("./grad_edit/grad_y.vtk", Y);
  handle.see_coord_grad_fields("./grad_edit/grad_z.vtk", Z);

  cout << "[info] done\n";
  return 0;
}
