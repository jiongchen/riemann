#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <jtflib/mesh/io.h>
#include <zjucad/matrix/io.h>

#include "src/spin_trans.h"

using namespace std;
using namespace zjucad::matrix;
using namespace riemann;
namespace po=boost::program_options;
using namespace Eigen;

namespace test_spin {
struct argument {
  string input_mesh;
  string curv_file;
  string output_folder;
};
}

static int read_curvature_change(const char *filename, matd_t &delta) {
  ifstream ifs(filename);
  if ( ifs.fail() ) {
    cerr << "[info] cannot open " << filename << endl;
    return __LINE__;
  }
  cout << "[info] load curvature change\n";
  size_t fid; double dv;
  while ( ifs >> fid >> dv ) {
    cout << "\tface: " << fid << " delta: " << dv << endl;
    delta[fid] = dv;
  }
  ifs.close();
  return 0;
}

int main(int argc, char *argv[])
{
  po::options_description desc("available options");
  desc.add_options()
      ("help,h", "produce help message")
      ("input_mesh,i", po::value<string>(), "input mesh")
      ("curv_file,c", po::value<string>(), "curvature change")
      ("output_folder,o", po::value<string>(), "output folder")
      ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if ( vm.count("help") ) {
    cout << desc << endl;
    return 1;
  }
  test_spin::argument args; {
    args.input_mesh = vm["input_mesh"].as<string>();
    args.curv_file = vm["curv_file"].as<string>();
    args.output_folder = vm["output_folder"].as<string>();
  }
  if ( !boost::filesystem::exists(args.output_folder) )
    boost::filesystem::create_directory(args.output_folder);

  // Input
  mati_t tris; matd_t nods;
  jtf::mesh::load_obj(args.input_mesh.c_str(), tris, nods);
  matd_t delta = zeros<double>(tris.size(2), 1);
  read_curvature_change(args.curv_file.c_str(), delta);

  // R^3->Im H
  matd_t nodx = zeros<double>(4, nods.size(2));
  nodx(colon(1, 3), colon()) = nods;
  matd_t xnew = nodx;

  spin_trans solver(tris, nodx);
  solver.set_curvature_change(delta);
  solver.deform(xnew);

  // Output
  char outfile[256];
  sprintf(outfile, "%s/origin.obj", args.output_folder.c_str());
  jtf::mesh::save_obj(outfile, tris, nods);
  nods = xnew(colon(1, 3), colon());
  sprintf(outfile, "%s/spin.obj", args.output_folder.c_str());
  jtf::mesh::save_obj(outfile, tris, nods);

  cout << "[info] done\n";
  return 0;
}
