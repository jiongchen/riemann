#include <iostream>
#include <fstream>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/util.h>
#include <zjucad/ptree/ptree.h>
#include <zjucad/matrix/io.h>

#include "src/polycube.h"
#include "src/vtk.h"
#include "src/blas_lapack.h"
#include <zjucad/matrix/lapack.h>

using namespace std;
using namespace riemann;
using namespace zjucad::matrix;

extern "C" {

  void surf_normal_align_(double *val, const double *x, const double *eps);
  void surf_normal_align_jac_(double *gra, const double *x, const double *eps);
  void surf_normal_align_hes_(double *hes, const double *x, const double *eps);
  
  void tet_distortion_(double *val, const double *x, const double *D, const double *R);
  void tet_distortion_jac_(double *jac, const double *x, const double *D, const double *R);
  void tet_distortion_hes_(double *hes, const double *x, const double *D, const double *R);

  void triangle_area_(double *val, const double *x);
  void triangle_area_jac_(double *jac, const double *x);

  void tri_area_normal_(double *out, const double *x);
  void tri_area_normal_jac_(double *out, const double *x);

}

namespace riemann {
  void area_normal_align_hes(const double *x, const double eps, const int id, double *out);
}

static void test_calc_tri_area() {
  matd_t nods = M_PI*rand(3, 3);
  double area = 0;
  triangle_area_(&area, &nods[0]);
  cout << "area diff: " << area-jtf::mesh::cal_face_area(nods) << endl;
}

static void test_calc_distortion() {
  matd_t nods = rand(3, 4);
  matd_t dm = nods(colon(), colon(1, 3))-nods(colon(), 0)*ones<double>(1, 3);
  inv(dm);

  matd_t F = (nods(colon(), colon(1, 3))-nods(colon(), 0)*ones<double>(1, 3))*dm;
  matd_t U, S, VT;
  svd(F, U, S, VT);
  matd_t R = U*VT;
  double value = 0;
  tet_distortion_(&value, &nods[0], &dm[0], &R[0]);
  cout << "distortion: " << value << endl;
}

static void test_calc_normal_align() {
  matd_t nods = sqrt(2)*eye<double>(3);
  nods(0, 0) = 0;
  double value = 0, eps = 0;
  surf_normal_align_(&value, &nods[0], &eps);
  cout << "normal align: " << value << endl;
}

static void test_lapack_eig() {
  const size_t dim = 5;
  
  matd_t A = rand(dim, dim);
  matd_t ATA = trans(A)+A;
  matd_t cache = ATA;
    
  matd_t e(dim, 1);
  eig(ATA, e);
  matd_t diag(dim, dim);
  for (size_t i = 0; i < dim; ++i)
    diag(i, i) = e[i];

  cout << "error: " << norm(cache-ATA*temp(diag*trans(ATA))) << endl;
}

static void test_component_wise_hes() {
  matd_t A = 1000*rand(3, 3);
  const double eps = 100;

  matd_t H(9, 9);
  surf_normal_align_hes_(&H[0], &A[0], &eps);
    
  matd_t Ht = zeros<double>(9, 9);
  for (int i = 0; i < 3; ++i) {
    matd_t tmp = zeros<double>(9, 9);
    area_normal_align_hes(&A[0], eps, i, &tmp[0]);
    Ht += tmp;
  }
  cout << "Hessian diff rel error: " << norm(H-Ht)/norm(H) << endl;  
}

static void test_area_normal() {
  cout << "---------------- TEST CASE ----------------" << endl;
  matd_t x = 10*rand(3, 3);
  matd_t n = zeros<double>(3, 1);
  const double eps = 100;

  tri_area_normal_(&n[0], &x[0]);
  double a1 = sqrt(n[0]*n[0]+eps)+sqrt(n[1]*n[1]+eps)+sqrt(n[2]*n[2]+eps);
  double a2 = 0;
  surf_normal_align_(&a2, &x[0], &eps);
  cout << "a1 value: " << a1 << endl << "a2 value: " << a2 << endl;
  cout << "align diff rel error: " << (a1-a2)/a2 << endl;

  matd_t g1 = zeros<double>(9, 1);
  matd_t jac = zeros<double>(3, 9);
  tri_area_normal_jac_(&jac[0], &x[0]);
  for (int i = 0; i < 3; ++i) {
    g1 += n[i]/sqrt(n[i]*n[i]+eps)*trans(jac(i, colon()));
  }
  matd_t g2 = zeros<double>(9, 1);
  surf_normal_align_jac_(&g2[0], &x[0], &eps);
  cout << "g1 norm: " << norm(g1) << endl << "g2 norm: " << norm(g2) << endl;
  cout << "align grad rel error: " << norm(g1-g2)/norm(g2) << endl;
}

int main(int argc, char *argv[])
{
#if 0
  srand(time(NULL));
  test_calc_tri_area();
  test_calc_distortion();
  test_calc_normal_align();
  test_lapack_eig();
  test_component_wise_hes();
  test_area_normal();
  return __LINE__;
#endif
  
  boost::property_tree::ptree pt;
  zjucad::read_cmdline(argc, argv, pt);

  const string outdir = pt.get<string>("outdir.value");
  
  mati_t tets; matd_t nods;
  jtf::mesh::tet_mesh_read_from_vtk(pt.get<string>("mesh.value").c_str(), &nods, &tets); {
    string outfile = outdir+string("/orig.vtk");
    ofstream ofs(outfile);
    tet2vtk(ofs, &nods[0], nods.size(2), &tets[0], tets.size(2));
    ofs.close();
  }

  shared_ptr<polycube_solver> polycube = make_shared<polycube_solver>(tets, nods, pt);

  cout << "[INFO] deform..." << endl;
  matd_t param_nods = nods;
  polycube->deform(param_nods);

  string outfile = outdir+string("/polycube.vtk");
  ofstream ofs(outfile);
  tet2vtk(ofs, &param_nods[0], param_nods.size(2), &tets[0], tets.size(2));
  ofs.close();
  
  cout << "[Info] done." << endl;
  return 0;
}
