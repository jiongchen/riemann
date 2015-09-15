#include <iostream>
#include <boost/filesystem.hpp>

#include "src/gradient_deform.h"

using namespace std;
using namespace riemann;

#define X 0
#define Y 1
#define Z 2

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
//  handle.see_coord_grad_fields("./grad_edit/grad_x.vtk", X);
//  handle.see_coord_grad_fields("./grad_edit/grad_y.vtk", Y);
//  handle.see_coord_grad_fields("./grad_edit/grad_z.vtk", Z);

//  handle.scale_grad_fields(1.5);
//  handle.reverse_grad_fields();
//  handle.see_coord_grad_fields("./grad_edit/scale_grad_x.vtk", X);
//  handle.see_coord_grad_fields("./grad_edit/scale_grad_y.vtk", Y);
//  handle.see_coord_grad_fields("./grad_edit/scale_grad_z.vtk", Z);

  vector<size_t> fixIDs{6976};
  handle.set_fixed_verts(fixIDs);
  vector<size_t> editIDs{1, 100, 200};
  handle.edit_boundary(editIDs);

  handle.propagate_transform();
  handle.see_harmonic_field("./grad_edit/hf.vtk");
  handle.deform();
  handle.save_deformed_model("./grad_edit/deform.obj");

  cout << "[info] done\n";
  return 0;
}
