#include <iostream>
#include <jtflib/mesh/io.h>
#include <boost/filesystem.hpp>

#include "src/arpaca.h"
#include "src/cotmatrix.h"

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
  SymmetricEigenSolver<double> sol = Solve(L, 20, ALGEBRAIC_SMALLEST);
  cout << sol.eigenvalues() << endl;

  cout << "[info] done\n";
  return 0;
}
