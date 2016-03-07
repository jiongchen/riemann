#include <iostream>
#include <unordered_set>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>

using namespace std;
using namespace zjucad::matrix;

int main(int argc, char *argv[])
{
  if ( argc != 3 ) {
    cerr << "#usage: " << argv[0] << " mesh.obj bnd_vert.txt\n";
    return __LINE__;
  }
  matrix<size_t> tris;
  matrix<double> nods;
  jtf::mesh::load_obj(argv[1], tris, nods);

  unordered_set<size_t> vert;
  shared_ptr<jtf::mesh::edge2cell_adjacent> e2c
      (jtf::mesh::edge2cell_adjacent::create(tris));
  for (auto &elem : e2c->edges_) {
    if ( e2c->is_boundary_edge(e2c->query(elem.first, elem.second)) ) {
      vert.insert(elem.first);
      vert.insert(elem.second);
    }
  }

  ofstream os(argv[2]);
  for (auto &vid : vert) {
    os << vid << endl;
  }
  os.close();
  cout << "[info] done\n";
  return 0;
}
