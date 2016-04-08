#include <iostream>
#include <boost/filesystem.hpp>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>

#include "src/json.h"
#include "igl/principal_curvature.h"
#include "igl/readOBJ.h"
#include "src/vtk.h"
#include "src/write_vtk.h"
#include "src/diffuse_dihedral_rot.h"
#include "src/dual_graph.h"

using namespace std;
using namespace Eigen;
using namespace riemann;
using namespace zjucad::matrix;

void tri2tet(const mati_t &tris, const matd_t &v_tri, mati_t &tets, matd_t &v_tet) {
  tets.resize(4, tris.size(2));
  v_tet.resize(3, v_tri.size(2)+tris.size(2));
  tets(colon(0, 2), colon()) = tris;
  tets(3, colon()) = colon(v_tri.size(2), v_tet.size(2)-1);
  v_tet(colon(), colon(0, v_tri.size(2)-1)) = v_tri;
#pragma omp parallel for
  for (size_t i = 0; i < tris.size(2); ++i) {
    matd_t vert = v_tri(colon(), tris(colon(), i));
    matd_t n = cross(vert(colon(), 1)-vert(colon(), 0), vert(colon(), 2)-vert(colon(), 0));
    v_tet(colon(), i+v_tri.size(2)) = vert(colon(), 0)+n/norm(n);
  }
}

extern "C" {
void calc_dihedral_angle_(double *value, const double *x);
}

void simple_test() {
  {
  Matrix<double, 3, 4> x;
  x << 1, 0, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, -1;
  cout << x << endl;
  double val = 1.0;
  calc_dihedral_angle_(&val, x.data());
  cout << val << endl;
  }

  {Matrix<double, 3, 4> x;
    x << 1, 0, 0, -1,
        0, 0, 1, 0,
        0, 0, 0, 0;
    cout << x << endl;
    double val = 1.0;
    calc_dihedral_angle_(&val, x.data());
    cout << val << endl;
    }
}

int main(int argc, char *argv[])
{
//  simple_test();
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

  shared_ptr<edge2cell_adjacent> ec;
  shared_ptr<Graph> g;
  build_tri_mesh_dual_graph(tris, ec, g);
  tree_t mst;
  get_minimum_spanning_tree(g, mst);

  const size_t root_face = json["root_face"].asUInt();
  vector<Matrix3d> rotation(tris.size(2));
  diffuse_rotation(tris, nods_prev, nods_curr, root_face, mst, rotation);

//  cout << "interrupt\n";
//  return 0;

  mati_t tets; matd_t tetv_prev, tetv_curr;
  tri2tet(tris, nods_prev, tets, tetv_prev);
  tri2tet(tris, nods_curr, tets, tetv_curr);
  {
    ofstream ofs(json["outdir"].asString()+"/tet_rest.vtk");
    tet2vtk(ofs, &tetv_prev[0], tetv_prev.size(2), &tets[0], tets.size(2));
    ofs.close();
  }
  {
    ofstream ofs(json["outdir"].asString()+"/tet_curr.vtk");
    tet2vtk(ofs, &tetv_curr[0], tetv_curr.size(2), &tets[0], tets.size(2));
    ofs.close();
  }

  diffuse_arap_solver solver(tets, tetv_prev, rotation);
  for (int i = 0; i < json["handles"].size(); ++i) {
    const size_t id = json["handles"][i]["id"].asUInt();
    solver.pin_down_vert(id, &tetv_curr(0, id), &tetv_curr[0]);
  }
  solver.solve(&tetv_curr[0]);

  {
    ofstream ofs(json["outdir"].asString()+"/tet_recover.vtk");
    tet2vtk(ofs, &tetv_curr[0], tetv_curr.size(2), &tets[0], tets.size(2));
    ofs.close();
  }
  {
    ofstream ofs(json["outdir"].asString()+"/tri_recover.vtk");
    tri2vtk(ofs, &tetv_curr[0], nods_curr.size(2), &tris[0], tris.size(2));
    ofs.close();
  }

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
