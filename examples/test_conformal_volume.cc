#include <iostream>
#include <fstream>
#include <jtflib/mesh/io.h>
#include <boost/filesystem.hpp>

#include "src/json.h"
#include "src/conformal_volume.h"
#include "src/vtk.h"

using namespace std;
using namespace riemann;
using namespace zjucad::matrix;

int main(int argc, char *argv[])
{
  if ( argc != 2 ) {
    cerr << "# usage: test_conformal_volume conf.json\n";
    return __LINE__;
  }
  Json::Reader reader;
  Json::Value json;
  ifstream ifs(argv[1]);
  if ( ifs.fail() ) {
    cerr << "[Error] can't open " << argv[1] << endl;
    return __LINE__;
  }
  if ( !reader.parse(ifs, json) ) {
    cerr << "[Error] " << reader.getFormattedErrorMessages() << endl;
    return __LINE__;
  }
  ifs.close();
  boost::filesystem::create_directories(json["output_dir"].asString());

  mati_t tets; matd_t nods;
  jtf::mesh::tet_mesh_read_from_zjumat(json["mesh"].asString().c_str(), &nods, &tets);
  matd_t verts = zeros<double>(4, nods.size(2));
  verts(colon(1, 3), colon()) = nods;

  conformal_volume cv(tets, verts);
  // parse charge
  for (unsigned int i = 0; i < json["charges"].size(); ++i) {
    const double pos[3] = {json["charges"][i]["pos"][0].asDouble(),
                           json["charges"][i]["pos"][1].asDouble(),
                           json["charges"][i]["pos"][2].asDouble()};
    const double intensity = json["charges"][i]["intensity"].asDouble();
    printf("# charge %u: (%lf, %lf, %lf), %lf\n", i, pos[0], pos[1], pos[2], intensity);
    cv.set_charge(pos, intensity);
  }
  cv.solve_eigen_prob();

  char outfile[256];
  sprintf(outfile, "%s/orig.vtk", json["output_dir"].asString().c_str());
  ofstream os(outfile);
  tet2vtk(os, &nods[0], nods.size(2), &tets[0], tets.size(2));
  point_data(os, cv.u_.data(), cv.u_.size(), "charge");
  os.close();

//  sprintf(outfile, "%s/grad.vtk", json["output_dir"].asString().c_str());
//  cv.draw_gradient(outfile);

  cout << "[Info] done\n";
  return 0;
}
