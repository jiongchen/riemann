#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <jtflib/mesh/io.h>
#include <unordered_set>

#include "src/advanced_mips.h"
#include "src/vtk.h"

using namespace std;
using namespace riemann;
using namespace zjucad::matrix;
namespace po=boost::program_options;

namespace test_amips {
struct argument {
  string src_mesh;
  string ini_mesh;
  string pos_cons;
  string out_folder;
};
}

static int read_fix_vert(const char *filename, unordered_set<size_t> &fix_vert) {
  ifstream is(filename);
  if ( is.fail() ) {
    cerr << "[info] can not open " << filename << endl;
    return __LINE__;
  }
  size_t vid;
  fix_vert.clear();
  while ( is >> vid )
    fix_vert.insert(vid);
  cout << "[info] number of fixed vertices: " << fix_vert.size() << endl;
  return 0;
}

int main(int argc, char *argv[])
{
  po::options_description desc("Available options");
  desc.add_options()
      ("help,h", "produce help message")
      ("source_mesh,s", po::value<string>(), "source mesh file")
      ("initial_mesh,i", po::value<string>(), "initial mesh file")
      ("pos_cons,c", po::value<string>(), "constraint file")
      ("output_folder,o", po::value<string>(), "output folder")
      ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if ( vm.count("help") ) {
    cout << desc << endl;
    return __LINE__;
  }
  test_amips::argument args; {
    args.src_mesh   = vm["source_mesh"].as<string>();
    args.ini_mesh   = vm["initial_mesh"].as<string>();
    args.pos_cons   = vm["pos_cons"].as<string>();
    args.out_folder = vm["output_folder"].as<string>();
  }
  if ( !boost::filesystem::exists(args.out_folder) )
    boost::filesystem::create_directory(args.out_folder);

  char filename[256];
  mati_t xz(2, 1); xz[0] = 0; xz[1] = 2; // y is zero by default
  mati_t tris; matd_t nods, nods0; {
    mati_t _tris; matd_t _nods, _nods0;
    jtf::mesh::load_obj(args.src_mesh.c_str(), _tris, _nods);
    jtf::mesh::load_obj(args.ini_mesh.c_str(), _tris, _nods0);
    sprintf(filename, "%s/source.obj", args.out_folder.c_str());
    jtf::mesh::save_obj(filename, _tris, _nods);
    sprintf(filename, "%s/initial.obj", args.out_folder.c_str());
    jtf::mesh::save_obj(filename, _tris, _nods0);
    tris = _tris;
    nods = _nods(xz, colon());
    nods0 = _nods0(xz, colon());
    cout << "[info] face number: " << tris.size(2) << endl;
    cout << "[info] vert number: " << _nods.size(2) << endl;
  }
  unordered_set<size_t> fix_vert;
  read_fix_vert(args.pos_cons.c_str(), fix_vert);

  mips_deformer_2d solver(tris, nods);
  solver.set_fixed_vert(fix_vert);
//  solver.unit_test();

  matd_t dir(2, 1);
  dir[0] = 0;
  dir[1] = 1;
//  for (size_t k = 0; k < 10; ++k) {
//    nods0(colon(), 22) += 0.05*dir;
//    solver.deform(&nods0[0], 1000);
//  }
  solver.deform(&nods0[0], 1001);

  sprintf(filename, "%s/deform.vtk", args.out_folder.c_str());
  ofstream os(filename); {
    matd_t _nods0 = zeros<double>(3, nods.size(2));
    _nods0(xz, colon()) = nods0;
    tri2vtk(os, &_nods0[0], _nods0.size(2), &tris[0], tris.size(2));
  }
  cout << "[info] done\n";
  return 0;
}
