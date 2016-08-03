#include <iostream>
#include <fstream>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <zjucad/ptree/ptree.h>

#include "src/polycube.h"
#include "src/vtk.h"

using namespace std;
using namespace riemann;
using namespace zjucad::matrix;

int main(int argc, char *argv[])
{
  boost::property_tree::ptree pt;
  zjucad::read_cmdline(argc, argv, pt);

  const string outdir = pt.get<string>("outdir.value");
  
  mati_t tets; matd_t nods;
  jtf::mesh::tet_mesh_read_from_vtk(pt.get<string>("mesh.value").c_str(), &nods, &tets); {
    string outfile = outdir+string("/orig.vtk");
    ofstream ofs(outfile);
    tet2vtk(ofs, &nods[0], nods.size(2), &tets[0], tets.size(2));
    ofs.close();
  }

  shared_ptr<polycube_solver> polycube = make_shared<polycube_solver>(tets, nods, pt);
#if 0
  polycube->unit_test();
#endif

  matd_t param_nods = nods;
  polycube->deform(param_nods);

  string outfile = outdir+string("/polycube.vtk");
  ofstream ofs(outfile);
  tet2vtk(ofs, &param_nods[0], param_nods.size(2), &tets[0], tets.size(2));
  ofs.close();
  
  cout << "[Info] done." << endl;
  return 0;
}
