#include <iostream>
#include <jtflib/mesh/io.h>
#include <boost/filesystem.hpp>

#include "src/arpaca.h"
#include "src/cotmatrix.h"
#include "src/vtk.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;
using namespace arpaca;

int main(int argc, char *argv[])
{
  if ( argc != 2 ) {
    cerr << "# Usage: " << argv[0] << " model.obj\n";
    return __LINE__;
  }
  boost::filesystem::create_directory("./lap_eigen");

  matrix<size_t> tris;
  matrix<double> nods;
  jtf::mesh::load_obj(argv[1], tris, nods);
  SparseMatrix<double> L;
  surfparam::cotmatrix(tris, nods, 1, &L);
  L = -L;
  SymmetricEigenSolver<double> sol = arpaca::Solve(L, 50, ALGEBRAIC_SMALLEST);
  cout << sol.eigenvalues() << endl;

  for (size_t i = 0; i < 50; ++i) {
    char out[256];
    sprintf(out, "./lap_eigen/mesh_wave_%zu.vtk", i);
    ofstream os(out);
    tri2vtk(os, &nods[0], nods.size(2), &tris[0], tris.size(2));
    VectorXd wave = sol.eigenvectors().col(i);
    point_data(os, &wave[0], wave.rows(), "eigen_func", "eigen_func");
    os.close();
  }

  cout << "[info] done\n";
  return 0;
}
