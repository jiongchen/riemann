#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <jtflib/mesh/io.h>

#include "src/local_edit.h"
#include "src/vtk.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;
using namespace riemann;

int main(int argc, char *argv[])
{
  if ( argc != 2 ) {
    cerr << "# usage: prog config.json\n";
    return __LINE__;
  }
  ifstream ifs(argv[1]);
  if ( ifs.fail() ) {
    cerr << "[ERROR] fail to open " << argv[1] << endl;
    return __LINE__;
  }
  Json::Value json;
  Json::Reader reader;
  bool flag = reader.parse(ifs, json);
  if ( !flag ) {
    cerr << "[ERROR] " << reader.getFormattedErrorMessages() << "\n";
    return __LINE__;
  }
  ifs.close();

  char outfile[256];
  boost::filesystem::create_directories(json["outdir"].asString());

  mati_t quad; matd_t nods, delta;
  jtf::mesh::load_obj(json["mesh"].asString().c_str(), quad, nods);
  sprintf(outfile, "%s/rest_pose.obj", json["outdir"].asString().c_str());
  jtf::mesh::save_obj(outfile, quad, nods);
  delta = zeros<double>(nods.size(1), nods.size(2));

  // editing and solve
  constrained_mesh_editor editor(quad, nods);
  editor.set_handles(json["handles"]);

  editor.deform(&delta[0]);

  matd_t nods_new = nods+delta;
  ofstream os(json["outdir"].asString()+string("/deform.vtk"));
  quad2vtk(os, &nods_new[0], nods_new.size(2), &quad[0], quad.size(2));

  cout << "[INFO] done\n";
  return 0;
}
