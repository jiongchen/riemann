#include <iostream>
#include <boost/filesystem.hpp>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>

#include "src/json.h"
#include "igl/principal_curvature.h"
#include "igl/readOBJ.h"
#include "src/vtk.h"
#include "src/write_vtk.h"

using namespace std;
using namespace Eigen;
using namespace riemann;
using namespace zjucad::matrix;
using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;
using jtf::mesh::edge2cell_adjacent;

extern "C" {
void calc_dihedral_angle_(double *val, const double *x);
void calc_dihedral_angle_jac_(double *jac, const double *x);
}

static void get_edge_elem(const mati_t &tris, mati_t &edge) {
  edge2cell_adjacent *e2c = edge2cell_adjacent::create(tris, false);
  edge.resize(2, e2c->edges_.size());
#pragma omp parallel
  for (size_t i = 0; i < edge.size(2); ++i) {
    edge(0, i) = e2c->edges_[i].first;
    edge(1, i) = e2c->edges_[i].second;
  }
  delete e2c;
}

static void get_diam_elem(const mati_t &tris, mati_t &diam) {
  edge2cell_adjacent *ea = edge2cell_adjacent::create(tris, false);
  mati_t bd_ed_id;
  get_boundary_edge_idx(*ea, bd_ed_id);
  diam.resize(4, ea->edges_.size()-bd_ed_id.size());
  for(size_t ei = 0, di = 0; ei < ea->edges_.size(); ++ei) {
    pair<size_t, size_t> nb_tr_id = ea->edge2cell_[ei];
    if( ea->is_boundary_edge(nb_tr_id) ) continue;
    diam(colon(1, 2), di) = ea->get_edge(ei);
    // orient
    bool need_swap = true;
    for(size_t k = 0; k < 3; ++k) {
      if( diam(1, di) == tris(k, nb_tr_id.first) ) {
        if( diam(2, di) != tris((k+1)%3, nb_tr_id.first) )
          need_swap = false;
      }
    }
    if( need_swap )
      swap(diam(1, di), diam(2, di));
    diam(0, di) = zjucad::matrix::sum(tris(colon(), nb_tr_id.first))
        - zjucad::matrix::sum(diam(colon(1, 2), di));
    diam(3, di) = zjucad::matrix::sum(tris(colon(), nb_tr_id.second))
        - zjucad::matrix::sum(diam(colon(1, 2), di));
    ++di;
  }
  delete ea;
}

int main(int argc, char *argv[])
{
  if ( argc != 2 ) {
    cerr << "#usage: ./test_ani_compression config.json\n";
    return __LINE__;
  }
  Json::Reader reader;
  Json::Value json;
  ifstream ifs(argv[1]);
  if ( ifs.fail() ) {
    cerr << "[Error] can not open " << argv[1] << endl;
    return __LINE__;
  }
  if ( !reader.parse(ifs, json) ) {
    cerr << "[Error] " << reader.getFormattedErrorMessages() << endl;
    return __LINE__;
  }
  ifs.close();

  mati_t tris; matd_t nods_prev, nods_curr;
  jtf::mesh::load_obj(json["mesh_prev"].asString().c_str(), tris, nods_prev);
  jtf::mesh::load_obj(json["mesh_curr"].asString().c_str(), tris, nods_curr);


  cout << "[Info] done\n";
  return 0;
}

//int main(int argc, char *argv[])
//{
//  if ( argc != 2 ) {
//    cerr << "#usage: ./test_ani_compression config.json\n";
//    return __LINE__;
//  }
//  Json::Reader reader;
//  Json::Value json;
//  ifstream ifs(argv[1]);
//  if ( ifs.fail() ) {
//    cerr << "[Error] can not open " << argv[1] << endl;
//    return __LINE__;
//  }
//  if ( !reader.parse(ifs, json) ) {
//    cerr << "[Error] " << reader.getFormattedErrorMessages() << endl;
//    return __LINE__;
//  }
//  ifs.close();

//  Matrix<int, -1, -1, RowMajor> tris;
//  Matrix<double, -1, -1, RowMajor> nods;
//  igl::readOBJ(json["mesh"].asString(), nods, tris);
//  Matrix<double, -1, -1, RowMajor> pd1, pd2;
//  VectorXd pv1, pv2;
//  igl::principal_curvature(nods, tris, pd1, pd2, pv1, pv2);

//  double scale = json["scale"].asDouble();
//  Matrix<double, -1, -1, RowMajor> cd1 = scale*pv1.asDiagonal()*pd1;
//  Matrix<double, -1, -1, RowMajor> cd2 = scale*pv2.asDiagonakl()*pd2;

//  for (int i = 0; i < pv1.size(); ++i)
//    printf("node: %d, min: %lf, max: %lf\n", i, pv1(i), pv2(i));

//  string outdir = json["outdir"].asString();
//  boost::filesystem::create_directories(outdir);
//  char buffer[256];
//  sprintf(buffer, "%s/model.vtk", outdir.c_str());
//  ofstream ofs(buffer);
//  tri2vtk(ofs, nods.data(), nods.rows(), tris.data(), tris.rows());
//  ofs.close();
//  sprintf(buffer, "%s/lam_min.vtk", outdir.c_str());
//  draw_vert_direct_field(buffer, nods.data(), nods.rows(), cd1.data());
//  sprintf(buffer, "%s/lam_max.vtk", outdir.c_str());
//  draw_vert_direct_field(buffer, nods.data(), nods.rows(), cd2.data());

//  cout << "[Info] done\n";
//  return 0;
//}
