#include <iostream>
#include <boost/filesystem.hpp>

#include "src/wave_constructor.h"

using namespace std;
using namespace riemann;

int main(int argc, char *argv[])
{
  if ( argc != 2 ) {
    cerr << "# usage: " << argv[0] << " model.obj\n";
    return __LINE__;
  }
  boost::filesystem::create_directory("./standing_wave");

  wave_constructor wc;
  wc.load_model_from_obj(argv[1]);
  wc.init();
  wc.solve_phase_transition();
  wc.prepare();
//  wc.solve_wave_value();

  wc.save_wave_to_vtk("./standing_wave/wave.vtk");

  cout << "[info] done\n";
  return 0;
}
