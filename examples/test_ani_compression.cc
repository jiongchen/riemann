#include <iostream>
#include <boost/filesystem.hpp>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>

#include "src/json.h"
#include "src/vtk.h"
#include "src/diffuse_dihedral_rot.h"
#include "src/dual_graph.h"
#include "src/config.h"
//#include "igl/principal_curvature.h"
//#include "igl/readOBJ.h"
//#include "src/write_vtk.h"

using namespace std;
using namespace Eigen;
using namespace riemann;
using namespace zjucad::matrix;

typedef signed char byte; // [-128, 127]

static void quantize(const vector<double> &dat, const pair<double, double> range,
                     const pair<byte, byte> bound, vector<byte> &qdat) {
  if ( qdat.size() != dat.size() )
    qdat.resize(dat.size());
  for (size_t i = 0; i < qdat.size(); ++i) {
    qdat[i] = floor(bound.first+(bound.second-bound.first)/(range.second-range.first)*(dat[i]-range.first)+0.5);
  }
}

static void dequantize(const vector<byte> &qdat, const pair<double, double> range,
                       const pair<byte, byte> bound, vector<double> &dat) {
  if ( dat.size() != qdat.size() )
    dat.resize(qdat.size());
  for (size_t i = 0; i < dat.size(); ++i) {
    dat[i] = range.first+(range.second-range.first)/(bound.second-bound.first)*(qdat[i]-bound.first);
  }
}

static int write_quant_res_bin(const char *file, const vector<byte> &qdat) {
  ofstream ofs(file, ios::binary);
  if ( ofs.fail() ) {
    cerr << "[Error] can not open " << file << endl;
    return __LINE__;
  }
  for (size_t i = 0; i < qdat.size(); ++i)
    ofs.write((char *)&qdat[i], sizeof(byte));
  ofs.close();
  return 0;
}

static int read_quant_res_bin(const char *file, vector<byte> &qdat) {
  ifstream ifs(file, ios::binary);
  if ( ifs.fail() ) {
    cerr << "[Error] can not open " << file << endl;
    return __LINE__;
  }
  if ( !qdat.empty() )
    qdat.clear();
  byte buffer;
  while ( ifs.read((char *)&buffer, sizeof(byte)) )
    qdat.push_back(buffer);
  ifs.close();
  return 0;
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

  // INPUT
  mati_t tris; matd_t nods, nods_prev, nods_curr;
  jtf::mesh::load_obj(json["mesh_rest"].asString().c_str(), tris, nods);
  jtf::mesh::load_obj(json["mesh_prev"].asString().c_str(), tris, nods_prev);
  jtf::mesh::load_obj(json["mesh_curr"].asString().c_str(), tris, nods_curr);

  string outdir = json["outdir"].asString();
  boost::filesystem::create_directories(outdir);

  // BUILD SPANNING TREE OF DUAL GRAPH
  shared_ptr<edge2cell_adjacent> ec;
  shared_ptr<Graph> g;
  build_tri_mesh_dual_graph(tris, ec, g);
  tree_t mst;
  get_minimum_spanning_tree(g, mst);
  const size_t root_face = json["root_face"].asUInt();

  // ENCODE
  diffuse_arap_encoder encoder;
  vector<double> da;
  encoder.calc_delta_angle(tris, nods_prev, nods_curr, mst, root_face, da);
  const size_t num_data = da.size();
  const double min_da = *std::min_element(da.begin(), da.end()),
      max_da = *std::max_element(da.begin(), da.end());
  printf("[Info] max delta: %lf\n[Info] min delta: %lf\n", min_da, max_da);

#define QUANTIZE_AND_COMPRESS
#ifdef QUANTIZE_AND_COMPRESS
  // QUANTIZE: UNIFORM SCALE AND TRUNCATE
  const byte BITS = json["bits"].asInt();
  ASSERT(BITS <= 8);
  const byte bound = (2<<(BITS-2))-1;
  printf("[Info] quantization bound: [%d, %d]\n", -bound, bound);
  vector<byte> cda;
  quantize(da, make_pair(min_da, max_da), make_pair(-bound, bound), cda);

  // WRITE FOR ZIP
  string quan_bin = outdir+string("/quant.dat");
  write_quant_res_bin(quan_bin.c_str(), cda);

  // READ FROM UNZIP
  cda.clear();
  read_quant_res_bin(quan_bin.c_str(), cda);
  ASSERT(cda.size() == num_data);
  cout << "[Info] read quantized data size: " << cda.size() << endl;

  // DEQUANTIZE
  da.clear();
  dequantize(cda, make_pair(min_da, max_da), make_pair(-bound, bound), da);
#endif

  // DECODE
  diffuse_arap_decoder decoder(tris, nods);
  matd_t root_curr = nods_curr(colon(), tris(colon(), root_face));
  decoder.estimate_rotation(nods_prev, mst, root_face, root_curr, da);
  for (size_t i = 0; i < 3; ++i) {
    const size_t id = tris(i, root_face);
    decoder.pin_down_vert(id, &nods_curr(0, id));
  }
  matd_t rec_curr(3, nods.size(2));
  decoder.solve(rec_curr);

  // OUTPUT
  {
    ofstream ofs(json["outdir"].asString()+"/tri_curr.vtk");
    tri2vtk(ofs, &nods_curr[0], nods_curr.size(2), &tris[0], tris.size(2));
    ofs.close();
  }
  {
    ofstream ofs(json["outdir"].asString()+"/tri_recover.vtk");
    tri2vtk(ofs, &rec_curr[0], rec_curr.size(2), &tris[0], tris.size(2));
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
