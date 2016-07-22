#include <iostream>
#include <fstream>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <boost/program_options.hpp>

#include "src/def.h"
#include "src/volume_frame.h"
#include "src/vtk.h"
#include "src/write_vtk.h"

using namespace std;
using namespace riemann;
using namespace Eigen;
using namespace zjucad::matrix;
namespace po=boost::program_options;

static int write_tet_zyz(const char *filename, const double *zyz, const size_t elem_num) {
  ofstream ofs(filename);
  if ( ofs.fail() ) {
    cerr << "[Error] can not write to " << filename << endl;
    return __LINE__;
  }
  ofs << "3 " << elem_num << endl;
  for (size_t i = 0; i < elem_num; ++i)
    ofs << zyz[3*i+0] << " " << zyz[3*i+1] << " " << zyz[3*i+2] << endl;
  return 0;
}

static int write_tet_ff(const char *filename, const double *ff, const size_t elem_num) {
  ofstream ofs(filename);
  if ( ofs.fail() ) {
    cerr << "[Error] can not write to " << filename << endl;
    return __LINE__;
  }
  ofs << "9 " << elem_num << endl;
  for (size_t i = 0; i < elem_num; ++i)
    ofs << ff[9*i+0] << " " << ff[9*i+1] << " " << ff[9*i+2] << " "
        << ff[9*i+3] << " " << ff[9*i+4] << " " << ff[9*i+5] << " "
        << ff[9*i+6] << " " << ff[9*i+7] << " " << ff[9*i+8] << endl;
  return 0;
}

int main(int argc, char *argv[])
{
  po::options_description desc("Available options"); {
    desc.add_options()
        ("help,h", "produce help message")
        ("mesh,i",          po::value<string>(), "input mesh")
        ("sm_type",         po::value<string>(), "smooth type")
        ("abs_eps",         po::value<double>(), "abs approximation")
        ("ws",              po::value<double>(), "smooth weight")
        ("wo",              po::value<double>(), "orth weight")
        ("wp",              po::value<double>(), "boundary weight")
        ("epsf",            po::value<double>()->default_value(1e-8), "epsf")
        ("maxits",          po::value<size_t>(), "maximum iterations")
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
  jtf::mesh::tet_mesh_read_from_vtk(vm["mesh"].as<string>().c_str(), &nods, &tets);
  {    
    string out = vm["output_folder"].as<string>()+string("/tet.vtk");
    ofstream ofs(out);
    tet2vtk(ofs, &nods[0], nods.size(2), &tets[0], tets.size(2));
    ofs.close();
  }

  // GET INITAL VALUE
  const double ws = 0.01, wa = 1e3, epsf = 1e-8;
  const size_t maxits = 10;
  cross_frame_args args = {ws, wa, epsf, maxits};

  shared_ptr<cross_frame_opt> frame_opt = make_shared<cross_frame_opt>(tets, nods, args);

  cout << "[INFO] solve Laplacian\n";
  VectorXd Fs;
  frame_opt->solve_laplacian(Fs);

  cout << "[INFO] solve for initial zyz angles\n";
  VectorXd abc;
  frame_opt->solve_initial_frames(Fs, abc);

  cout << "[INFO] optimize frames\n";
  frame_opt->optimize_frames(abc);

  
  // OPTIMZE SMOOTHNESS
  smooth_args sm_args = {
    vm["abs_eps"].as<double>(),
    vm["ws"].as<double>(),
    vm["wo"].as<double>(),
    vm["wp"].as<double>(),
    vm["epsf"].as<double>(),
    vm["maxits"].as<size_t>()
  };

  shared_ptr<frame_smoother> smoother = make_shared<frame_smoother>(tets, nods, sm_args);

  cout << "[Info] Fix boundary and opt smooth energy\n";

  if ( vm["sm_type"].as<string>() == "SH" ) {
    cout << "\t**************** SH ****************\n";
    smoother->smoothSH(abc);

    string zyz_file = vm["output_folder"].as<string>()+string("/frames.txt");
    write_tet_zyz(zyz_file.c_str(), abc.data(), abc.size()/3);

  } else if ( vm["sm_type"].as<string>() == "L1" ) {
    cout << "\t**************** L1 ****************\n";
    VectorXd fmat;
    convert_zyz_to_mat(abc, fmat);
    smoother->smoothL1(fmat);

    string ff_file = vm["output_folder"].as<string>()+string("/frames.txt");
    write_tet_ff(ff_file.c_str(), fmat.data(), fmat.size()/9);

  } else {
    cerr << "[Error] unsupported smooth term.\n";
    return __LINE__;
  }
  
  cout << "[Info] done\n";
  return 0;
}
