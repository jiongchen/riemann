#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <jtflib/mesh/io.h>

#include "src/advanced_mips.h"
#include "src/vtk.h"

using namespace std;
using namespace riemann;
using namespace zjucad::matrix;
namespace po=boost::program_options;

struct argument {
  string src_mesh;
  string ini_mesh;
  string pos_cons;
  string out_folder;
};

int main(int argc, char *argv[])
{

  mati_t tris; matd_t nods; {
    mati_t _tris; matd_t _nods;
    jtf::mesh::load_obj(argv[1], _tris, _nods);
    tris = _tris;
    nods = _nods;
  }

  mips_deformer_2d solver;
  solver.unit_test();

  ofstream os(); {
    matd_t _nods = zeros<double>(3, nods.size(2));
//    tri2vtk(os, &_nods[0], _nods.size(2), &tris[0], tris.size(2));
  }
  cout << "[info] done\n";
  return 0;
}
