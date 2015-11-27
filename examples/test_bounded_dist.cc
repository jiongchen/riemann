#include <iostream>
#include <jtflib/mesh/io.h>
#include <boost/filesystem.hpp>
#include <zjucad/matrix/io.h>
#include <zjucad/matrix/itr_matrix.h>

#include "src/bounded_distortion.h"
#include "src/vtk.h"

using namespace std;
using namespace riemann;
using namespace Eigen;
using namespace zjucad::matrix;

int main(int argc, char *argv[])
{
  if ( argc != 3 ) {
    cerr << "# usage: prog rest.tet.vtk deform.tet.vtk\n";
    return __LINE__;
  }
  boost::filesystem::create_directory("./bds");

  mati_t tets; matd_t nods, nods0;
  jtf::mesh::tet_mesh_read_from_vtk(argv[1], &nods, &tets);
  jtf::mesh::tet_mesh_read_from_vtk(argv[2], &nods0, &tets);

  SparseMatrix<double> T;
  matd_t binv;
  calc_tet_base_inv(tets, nods, binv);
  calc_tet_defo_grad_map(tets, binv, &T);

  // 1
  matd_t defoGrad(9, tets.size(2));
  Map<VectorXd>(&defoGrad[0], defoGrad.size()) = T*Map<VectorXd>(&nods0[0], nods0.size());
  // 2
  matd_t defGrad1(9, tets.size(2));
  for (size_t i = 0; i < tets.size(2); ++i) {
    matd_t base = nods0(colon(), tets(colon(1, 3), i))-nods0(colon(), tets(0, i))*ones<double>(1, 3);
    itr_matrix<double *>(3, 3, &defGrad1(0, i)) = base*itr_matrix<double *>(3, 3, &binv(0, i));
  }
  cout << defoGrad(colon(), 0) << endl;
  cout << max(defoGrad-defGrad1) << endl;

  cout << "[info] done\n";
  return 0;
}
