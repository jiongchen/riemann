#include <iostream>
#include <fstream>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <boost/program_options.hpp>

#include "src/volume_frame.h"
#include "src/vtk.h"

using namespace std;
using namespace riemann;
using namespace Eigen;
namespace po=boost::program_options;

int main(int argc, char *argv[])
{
  po::options_description desc("Available options"); {
    desc.add_options()
        ("help,h", "produce help message")
        ("mesh,i", po::value<string>(), "input mesh")
        ("output_folder,o", po::value<string>(), "output folder")
        ;
  }
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if ( vm.count("help") ) {
    cout << desc << endl;
    return __LINE__;
  }

  mati_t tets; matd_t nods;
  jtf::mesh::tet_mesh_read_from_zjumat(vm["mesh"].as<string>().c_str(), &nods, &tets);

  shared_ptr<cross_frame_opt> frame_opt(cross_frame_opt::create(tets, nods));

  VectorXd Fs = VectorXd::Zero(9*nods.size(2));
  frame_opt->solve_smooth_sh_coeffs(Fs);

  VectorXd abc = VectorXd::Zero(3*nods.size(2));
  frame_opt->solve_initial_frames(abc);
  frame_opt->optimize_frames(abc);


  cout << "[Info] done\n";
  return 0;
}
