#include <iostream>
#include <boost/filesystem.hpp>

#include "src/wave_constructor.h"

using namespace std;
using namespace riemann;

int main(int argc, char *argv[])
{
  if ( argc != 4 ) {
    cerr << "# usage: " << argv[0] << " model.obj frame.ef feature.fl\n";
    return __LINE__;
  }
  boost::filesystem::create_directory("./standing_wave");

  wave_constructor wc;
  wc.load_model_from_obj(argv[1]);
  wc.load_frame_field(argv[2]);
  wc.load_feature_line(argv[3]);

  wc.vis_edge_frame_field("./standing_wave/X.vtk", "./standing_wave/Y.vtk", 0.01);
  wc.save_feature_to_vtk("./standing_wave/feature_line.vtk");

  wc.scale_frame_field(0.02);
  wc.solve_phase_transition();
  wc.prepare();
  wc.give_an_initial_value(1);
  wc.solve_wave();

  wc.save_wave_to_vtk("./standing_wave/wave.vtk");

  cout << "[info] done\n";
  return 0;
}
