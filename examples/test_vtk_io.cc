#include <iostream>
#include <boost/filesystem.hpp>
#include <jtflib/mesh/io.h>
#include <unordered_set>

#include "src/cotmatrix.h"
#include "src/write_vtk.h"
#include "src/util.h"
#include "src/grad_operator.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;
using namespace riemann;

int main(int argc, char *argv[])
{
  if ( argc != 2 ) {
    cerr << "# Uasage: " << argv[0] << " model.obj\n";
    return __LINE__;
  }
  boost::filesystem::create_directory("./vtk_io");

  matrix<size_t> tris;
  matrix<double> nods;
  jtf::mesh::load_obj(argv[1], tris, nods);

  // compute cotan weight laplacian matrix
  SparseMatrix<double> L;
  cotmatrix(tris, nods, 1, &L);
  cout << "[info] laplacian size: " << L.rows() << endl;

  // init a scalar field with some vertices evaluated 1.0
  VectorXd sf = VectorXd::Zero(nods.size(2));
  srand(time(NULL));
  vector<size_t> idx(2);
  for (size_t i = 0; i < idx.size(); ++i) {
    idx[i] = rand() % nods.size(2);
    sf[idx[i]] = 1.0-i%2;
  }
  draw_vert_value_to_vtk("./vtk_io/sf.vtk", nods.begin(), nods.size(2),
                         tris.begin(), tris.size(2), sf.data());

  // solve for a harmonic field
  unordered_set<size_t> fixDOF;
  for (auto &id : idx)
    fixDOF.insert(id);
  vector<size_t> g2l;
  build_global_local_mapping<size_t>(nods.size(2), fixDOF, g2l);
  VectorXd rhs = -L*sf;
  rm_spmat_col_row(L, g2l);
  rm_vector_row(rhs, g2l);
  SimplicialCholesky<SparseMatrix<double>> sol;
  sol.compute(L);
  ASSERT(sol.info() == Success);
  VectorXd df = sol.solve(rhs);
  ASSERT(sol.info() == Success);
  VectorXd Df = VectorXd::Zero(nods.size(2));
  rc_vector_row(df, g2l, Df);
  sf += Df;
  draw_vert_value_to_vtk("./vtk_io/harmonic_sf.vtk", nods.begin(), nods.size(2),
                         tris.begin(), tris.size(2), sf.data());

  // compute gradient
  SparseMatrix<double> G;
  calc_grad_operator(tris, nods, &G);
  VectorXd grad = G*sf;
  grad *= 0.01;
  draw_face_direct_field("./vtk_io/grad.vtk", nods.begin(), nods.size(2),
                         tris.begin(), tris.size(2), grad.data());

  cout << "[info] done\n";
  return 0;
}
