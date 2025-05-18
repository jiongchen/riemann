#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>

#include "src/frame_field_deform.h"

using namespace std;
using namespace riemann;

int main(int argc, char *argv[])
{
  if ( argc != 3 ) {
    cerr << "# usage: " << argv[0] << " model.obj cons.ff\n";
    return __LINE__;
  }
  boost::filesystem::create_directories("./result/ff_deform");

  boost::property_tree::ptree pt;
  pt.put("max_iter", 2000);
  pt.put("tolerance", 1e-12);
  pt.put("lambda", 0.1);
  pt.put("perturb", 0.1);

  frame_field_deform deformer(pt);
  deformer.load_mesh(argv[1]);
  deformer.save_original_mesh("./result/ff_deform/origin.obj");
  deformer.visualize_local_bases("./result/ff_deform/local_frame.vtk", 0.1);

  deformer.load_constraints(argv[2]);
  deformer.visualize_init_frames("./result/ff_deform/init_ff.vtk");
  deformer.visualize_frame_fields("./result/ff_deform/interp_ff.vtk");

  deformer.interp_frame_fields();
  bool spd = deformer.check_spd_tensor_fields();
  cout << "[INFO] spd tensor: " << spd << endl;
  deformer.visualize_tensor_fields("./result/ff_deform/tensor.vtk");

  deformer.precompute();
  deformer.deform();

  deformer.save_original_mesh("./result/ff_deform/origin_post.obj");
  deformer.save_deformed_mesh("./result/ff_deform/deform.obj");

  deformer.calc_defo_grad_oper();
  deformer.gen_cross_field();
  deformer.visualize_cross_fields("./result/ff_deform/cross.vtk");

  cout << "[INFO] done\n";
  return 0;
}
