#include <iostream>
#include <boost/filesystem.hpp>

#include "src/vec_field_deform.h"

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[])
{
  boost::filesystem::create_directories("./result/vel_field_deform");

  riemann::vel_field_deform def;
  def.load_model("./dat/beam.obj");

  // // for beam twist
  // Vector3d center(-6, 0, 0), n(1, 0, 0);
  // const double ri = 12.05;
  // const double ro = 14;
  // def.twist_deform(center, ri, ro, n, 400);

  // for beam bend
  Vector3d center(0, 0, 0), n(1, 0, 0), axis(0, 0, -1);
  const double ri = -2;
  const double ro = 2;
  def.bend_deform(center, ri, ro, axis, n, M_PI/2);

  def.save_model("./result/vel_field_deform/deform.obj");
  cout << "done\n";
  return 0;
}
