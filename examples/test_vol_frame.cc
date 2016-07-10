#include <iostream>
#include <fstream>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <boost/program_options.hpp>

#include "src/volume_frame.h"
#include "src/vtk.h"
#include "src/lbfgs_solve.h"
#include "src/sh_zyz_convert.h"
#include "src/write_vtk.h"
#include "src/geometry_extend.h"

using namespace std;
using namespace riemann;
using namespace Eigen;
namespace po=boost::program_options;

class test_func : public Functional<double>
{
public:
  size_t Nx() const {
    return 2;
  }
  int Val(const double *x, double *val) const {
    *val += 100*pow(x[0]+3,4) + pow(x[1]-3,4);
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    gra[0] += 400*pow(x[0]+3,3);
    gra[1] += 4*pow(x[1]-3,3);
    return 0;
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    return __LINE__;
  }
};

static int simple_test_on_lbfgs()
{
  shared_ptr<Functional<double>> F = make_shared<test_func>();
  double x[2] = {0};
  
  const double epsf = 0, epsx = 0;
  const size_t maxits = 0;
  lbfgs_solve(F, x, 2, epsf, epsx, maxits);

  cout << x[0] << " " << x[1] << endl;
  return 0;
}

static int write_tet_zyz(const char *filename, const MatrixXd &tet_zyz) {
  ofstream ofs(filename);
  if ( ofs.fail() ) {
    cerr << "[Error] can not write to " << filename << endl;
    return __LINE__;
  }
  ofs << tet_zyz.rows() << " " << tet_zyz.cols() << endl;
  for (size_t i = 0; i < tet_zyz.cols(); ++i)
    ofs << tet_zyz(0, i) << " " << tet_zyz(1, i) << " " << tet_zyz(2, i) << endl;
  return 0;
}

static int zyz_vert_to_tet(const mati_t &tets, const VectorXd &abc, MatrixXd& zyz) {
  Map<const MatrixXd> ABC(abc.data(), 3, abc.size()/3);
  for (size_t i = 0; i < tets.size(2); ++i) {
    zyz.col(i) = (ABC.col(tets(0, i))+ABC.col(tets(1, i))+ABC.col(tets(2, i))+ABC.col(tets(3, i)))/4.0;
  }
  return 0;
}

static void zyz_sh_convert_test() {
  srand(time(NULL));
  
  Vector3d zyz = Vector3d::Random();
  cout << "orig zyz: " << zyz.transpose() << endl;

  Matrix<double, 9, 1> sh;
  zyz_to_sh(zyz.data(), sh.data());
  cout << "SH coeffs: " << sh.transpose() << endl;

  sh_to_zyz(sh.data(), zyz.data(), 10);
  cout << "recover zyz: " << zyz.transpose() << endl;
}

extern "C" {
  void cubic_sym_align_(double *val, const double *, const double *, const double *area);
}

static void alignment_test() {
  Vector3d zyz = Vector3d::Random();
  Matrix3d frm = RZ(zyz[2])*RX(-M_PI/2)*RZ(zyz[1])*RX(M_PI/2)*RZ(zyz[0]);
  Vector3d axis = frm.col(2);
  Vector3d nzyz = Vector3d(-atan2(axis[1], axis[0]), -acos(axis[2]), 0);
  double value = 0, area = 1;
  
  cubic_sym_align_(&value, zyz.data(), nzyz.data(), &area);
  cout << value << endl << endl;

  zyz = Vector3d::Ones()*M_PI; //setZero();
  VectorXd sh(9);
  zyz_to_sh(zyz.data(), sh.data());
  cout << sh << endl;
}

int main(int argc, char *argv[])
{
  alignment_test();
  return __LINE__;
  
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

  string tet_out = vm["output_folder"].as<string>()+string("/tet.vtk");
  ofstream ofs(tet_out);
  tet2vtk(ofs, &nods[0], nods.size(2), &tets[0], tets.size(2));
  ofs.close();

  shared_ptr<cross_frame_opt> frame_opt(cross_frame_opt::create(tets, nods));

  cout << "[INFO] solve Laplacian\n";
  VectorXd Fs = VectorXd::Zero(9*nods.size(2));
  frame_opt->solve_laplacian(Fs);
  
  cout << "[INFO] solve for initial zyz angles\n";
  VectorXd abc = VectorXd::Zero(3*nods.size(2));
  frame_opt->solve_initial_frames(Fs, abc);

  cout << "[INFO] optimize frames\n";
  frame_opt->optimize_frames(abc);

  MatrixXd frames(9, nods.size(2));
  for (size_t i = 0; i < nods.size(2); ++i) {
    Matrix3d R = RZ(abc[3*i+2])*RY(abc[3*i+1])*RZ(abc[3*i+0]);
    Map<Matrix3d>(&frames(0, i)) = R;
  }
  frames *= 0.05;
  
  string xfile = vm["output_folder"].as<string>()+string("/x.vtk");
  MatrixXd rx = frames.block(0, 0, 3, frames.cols());
  draw_vert_direct_field(xfile.c_str(), &nods[0], nods.size(2), rx.data());

  string yfile = vm["output_folder"].as<string>()+string("/y.vtk");
  MatrixXd ry = frames.block(3, 0, 3, frames.cols());
  draw_vert_direct_field(yfile.c_str(), &nods[0], nods.size(2), ry.data());

  string zfile = vm["output_folder"].as<string>()+string("/z.vtk");
  MatrixXd rz = frames.block(6, 0, 3, frames.cols());
  draw_vert_direct_field(zfile.c_str(), &nods[0], nods.size(2), rz.data());
  
  // // interpolate frames on tet
  // MatrixXd tet_zyz(3, tets.size(2));
  // zyz_vert_to_tet(tets, abc, tet_zyz);

  // // write zyz
  // string zyz_file = vm["output_folder"].as<string>()+string("/zyz.txt");
  // write_tet_zyz(zyz_file.c_str(), tet_zyz);
  
  cout << "[Info] done\n";
  return 0;
}
