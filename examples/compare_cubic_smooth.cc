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

int main(int argc, char *argv[])
{
  po::options_description desc("Available options"); {
    desc.add_options()
        ("help,h", "produce help message")
        ("mesh,i",          po::value<string>(), "input mesh")
        ("abs_eps",         po::value<double>(), "abs approximation")
        ("ws",              po::value<double>(), "smooth weight")
        ("wo",              po::value<double>(), "orth weight")
        ("wp",              po::value<double>(), "boundary weight")
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
  const double ws = 0.5, wa = 1e3, epsf = 1e-8;
  const size_t maxits = 50;
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
    1e-8,
    vm["maxits"].as<size_t>()
  };

  shared_ptr<frame_smoother> smoother = make_shared<frame_smoother>(tets, nods, sm_args);

  cout << "[Info] Fix boundary and opt smooth energy\n";

  // smoother->smoothSH(abc);

  VectorXd fmat;
  convert_zyz_to_mat(abc, fmat);
  smoother->smoothL1(fmat);

  
  // // write frame vectors
  // VectorXd fmat;
  // convert_zyz_to_mat(abc, fmat);

  MatrixXd frames = Map<MatrixXd>(fmat.data(), 9, fmat.size()/9);
  const double len_scale = 0.075;
  frames *= len_scale;

  matd_t center(3, tets.size(2));
  for (size_t i = 0; i < tets.size(2); ++i)
    center(colon(), i) = nods(colon(), tets(colon(), i))*ones<double>(4, 1)/4.0;
  
  string xfile = vm["output_folder"].as<string>()+string("/x.vtk");
  MatrixXd rx = frames.block(0, 0, 3, frames.cols());
  draw_vert_direct_field(xfile.c_str(), &center[0], center.size(2), rx.data());

  string yfile = vm["output_folder"].as<string>()+string("/y.vtk");
  MatrixXd ry = frames.block(3, 0, 3, frames.cols());
  draw_vert_direct_field(yfile.c_str(), &center[0], center.size(2), ry.data());

  string zfile = vm["output_folder"].as<string>()+string("/z.vtk");
  MatrixXd rz = frames.block(6, 0, 3, frames.cols());
  draw_vert_direct_field(zfile.c_str(), &center[0], center.size(2), rz.data());

  // // write zyz
  // string zyz_file = vm["output_folder"].as<string>()+string("/zyz.txt");
  // write_tet_zyz(zyz_file.c_str(), abc.data(), abc.size()/3);
  
  cout << "[Info] done\n";
  return 0;
}
