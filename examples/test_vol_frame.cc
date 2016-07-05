#include <iostream>
#include <fstream>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <boost/program_options.hpp>

#include "src/volume_frame.h"
#include "src/vtk.h"

using namespace std;
using namespace riemann;
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

  string outfile = vm["output_folder"].as<string>()+string("/orig.vtk");
  ofstream ofs(outfile);
  tet2vtk(ofs, &nods[0], nods.size(2), &tets[0], tets.size(2));
  ofs.close();
  
  cout << "[Info] done\n";
  return 0;
}
