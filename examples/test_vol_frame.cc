#include <iostream>
#include <fstream>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <zjucad/ptree/ptree.h>

#include "src/def.h"
#include "src/volume_frame.h"
#include "src/vtk.h"
#include "src/lbfgs_solve.h"
#include "src/sh_zyz_convert.h"
#include "src/write_vtk.h"
#include "src/geometry_extend.h"

using namespace std;
using namespace riemann;
using namespace Eigen;
using namespace zjucad::matrix;

extern "C" {
  void cubic_sym_align_(double*, const double*, const double*, const double*);
  void poly_smooth_tet_(double *val, const double *abc, const double *stiff);
  void cubic_sym_smooth_tet_(double *val, const double *abc, const double *stiff);
}

class test_func : public Functional<double>
{
public:
  size_t Nx() const {
    return 2;
  }
  int Val(const double *x, double *val) const {
    *val += 100*pow(x[0]+3, 4)+pow(x[1]-3, 4);
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    gra[0] += 400*pow(x[0]+3, 3);
    gra[1] += 4*pow(x[1]-3, 3);
    return 0;
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    return __LINE__;
  }
};

static void test_lbfgs_opt() {
  cout << "---------------TEST CASE----------------" << endl;
  shared_ptr<Functional<double>> F = make_shared<test_func>();
  double x[2] = {0};
  
  const double epsf = 0, epsx = 0;
  const size_t maxits = 0;
  lbfgs_solve(F, x, 2, epsf, epsx, maxits);

  cout << x[0] << " " << x[1] << endl << endl;
}

static void test_alignment_energy() {
  cout << "---------------TEST CASE----------------" << endl;
  srand(time(NULL));
  Vector3d zyz = Vector3d::Random();
  Matrix3d frm = RZ(zyz[2])*RX(-M_PI/2)*RZ(zyz[1])*RX(M_PI/2)*RZ(zyz[0]);

  for (size_t j = 0; j < 3; ++j) {
    Vector3d axis = frm.col(j);
    Vector3d nzyz = Vector3d(-atan2(axis[1], axis[0]), -acos(axis[2]), 0);
    double value = 0, area = 1;
    cubic_sym_align_(&value, zyz.data(), nzyz.data(), &area);
    cout << "aixs " << j << ": " << value << endl;
  }
  cout << endl;
}

static void test_sh_to_zyz() {
  cout << "---------------TEST CASE----------------" << endl;
  srand(time(NULL));
  
  Vector3d zyz = Vector3d::Random();
  cout << "orig zyz: " << zyz.transpose() << endl;

  Matrix<double, 9, 1> sh = Matrix<double, 9, 1>::Zero();
  zyz_to_sh(zyz.data(), sh.data());
  cout << "SH coeffs: " << sh.transpose() << endl;
  double res = 0;
  sh_residual_(&res, zyz.data(), sh.data());
  cout << "residual: " << res << endl;
  
  zyz.setZero();
  sh_to_zyz(sh.data(), zyz.data(), 1000);
  cout << "recover zyz: " << zyz.transpose() << endl;
  res = 0;
  sh_residual_(&res, zyz.data(), sh.data());
  cout << "reconstructed residual: " << res << endl;
  zyz_to_sh(zyz.data(), sh.data());
  cout << "SH coeffs: " << sh.transpose() << endl << endl;
}

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

static int verify_sh_scale() {
  cout << "---------------TEST CASE----------------" << endl;
  srand(time(NULL));
  Vector3d abc = Vector3d::Random();
  cout << "zyz: " << abc.transpose() << endl;
  const double stiff = 1.0;
  double v0 = 0, v1 = 0;
  poly_smooth_tet_(&v0, abc.data(), &stiff);
  cubic_sym_smooth_tet_(&v1, abc.data(), &stiff);
  v1 *= 4*M_PI/(15*15*7);
  cout << "SH: " << v1 << endl;
  cout << "Poly: " << v0 << endl;
  cout << "SH/poly: " << v1/v0 << endl;
  cout << 16*M_PI/315 << endl;
  return 0;
}

int main(int argc, char *argv[])
{
#if 0
  test_lbfgs_opt();
  test_alignment_energy();
  test_sh_to_zyz();
  verify_sh_scale();
  return __LINE__;
#endif

  boost::property_tree::ptree pt;
  zjucad::read_cmdline(argc, argv, pt);
  
  const string out_folder = pt.get<string>("out_dir.value");
  
  mati_t tets; matd_t nods;
  jtf::mesh::tet_mesh_read_from_vtk(pt.get<string>("mesh.value").c_str(), &nods, &tets);
  {    
    string outfile = out_folder+string("/tet.vtk");
    ofstream ofs(outfile);
    tet2vtk(ofs, &nods[0], nods.size(2), &tets[0], tets.size(2));
    ofs.close();
  }
  
  shared_ptr<cross_frame_opt> frame_opt = make_shared<cross_frame_opt>(tets, nods, pt);

  cout << "[INFO] solve Laplacian\n";
  VectorXd Fs;
  frame_opt->solve_laplacian(Fs);

  cout << "[INFO] solve for initial zyz angles\n";
  VectorXd abc;
  frame_opt->solve_initial_frames(Fs, abc);

  cout << "[INFO] optimize frames\n";
  frame_opt->optimize_frames(abc);

  // write frame vectors
  MatrixXd frames; {
    VectorXd fmat;
    convert_zyz_to_mat(abc, fmat);
    frames = Map<MatrixXd>(fmat.data(), 9, fmat.size()/9);
    const double len_scale = 0.075;
    frames *= len_scale;
  }
  
  matd_t bc(3, tets.size(2));
  for (size_t i = 0; i < tets.size(2); ++i)
    bc(colon(), i) = nods(colon(), tets(colon(), i))*ones<double>(4, 1)/4.0;
  
  string xfile = out_folder+string("/x.vtk");
  MatrixXd rx = frames.block(0, 0, 3, frames.cols());
  draw_vert_direct_field(xfile.c_str(), &bc[0], bc.size(2), rx.data());

  string yfile = out_folder+string("/y.vtk");
  MatrixXd ry = frames.block(3, 0, 3, frames.cols());
  draw_vert_direct_field(yfile.c_str(), &bc[0], bc.size(2), ry.data());

  string zfile = out_folder+string("/z.vtk");
  MatrixXd rz = frames.block(6, 0, 3, frames.cols());
  draw_vert_direct_field(zfile.c_str(), &bc[0], bc.size(2), rz.data());

  // write zyz
  string zyz_file = out_folder+string("/zyz.txt");
  write_tet_zyz(zyz_file.c_str(), abc.data(), abc.size()/3);
  
  cout << "[Info] done\n";
  return 0;
}
