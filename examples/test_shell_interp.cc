#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <jtflib/mesh/io.h>

#include "src/shell.h"

using namespace std;
using namespace zjucad::matrix;
using namespace riemann;
namespace po=boost::program_options;

struct argument {
  string input_mesh;
  string output_folder;
  shell_args sa;
};

int main(int argc, char *argv[])
{
  po::options_description desc("available options");
  desc.add_options()
      ("help,h", "produce help message")
      ("input_mesh,i", po::value<string>(), "input mesh")
      ("output_folder,o", po::value<string>(), "output folder")
      ("ws", po::value<double>()->default_value(1e3), "weight of stretch")
      ("wb", po::value<double>()->default_value(1e0), "weight of bending")
      ("wp", po::value<double>()->default_value(1e4), "weight of position")
      ("max_iter,n", po::value<size_t>()->default_value(10000), "max iter")
      ("tolerance,e", po::value<double>()->default_value(1e-12), "tolerance")
      ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if ( vm.count("help") ) {
    cout << desc << endl;
    return 1;
  }
  argument args;
  args.input_mesh = vm["input_mesh"].as<string>();
  args.output_folder = vm["output_folder"].as<string>();
  args.sa.ws = vm["ws"].as<double>();
  args.sa.wb = vm["wb"].as<double>();
  args.sa.wp = vm["wp"].as<double>();
  args.sa.max_iter = vm["max_iter"].as<size_t>();
  args.sa.tolerance = vm["tolerance"].as<double>();

  if ( !boost::filesystem::exists(args.output_folder) )
    boost::filesystem::create_directory(args.output_folder);

  matrix<size_t> tris;
  matrix<double> nods;
  jtf::mesh::load_obj(args.input_mesh.c_str(), tris, nods);

  shell_deformer defo(tris, nods, args.sa); {
    size_t idx = 3;
    matd_t u = nods(colon(), idx);
    defo.fix_vert(idx, &u[0]);
  } {
    size_t idx = 0;
    matd_t u = nods(colon(), idx)-0.1*ones<double>(3, 1);
    u[1] += 0.4;
    defo.fix_vert(idx, &u[0]);
  }
  defo.prepare();
  defo.solve(&nods[0]);

  char outfile[256];
  sprintf(outfile, "%s/deform.obj", args.output_folder.c_str());
  jtf::mesh::save_obj(outfile, tris, nods);

  cout << "done\n";
  return 0;
}
