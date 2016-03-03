#include <iostream>
#include <boost/filesystem.hpp>

#include "igl/readOBJ.h"
#include "igl/writeOBJ.h"
#include "igl/boundary_loop.h"

#include "src/json.h"
#include "src/vtk.h"
#include "src/hole_filling.h"
#include "src/vtk.h"
#include "src/write_vtk.h"

using namespace std;
using namespace Eigen;
using namespace riemann;

int main(int argc, char *argv[])
{
  if ( argc != 2 ) {
    cerr << "# usage: prog conf.json\n";
    return __LINE__;
  }
  Json::Reader reader;
  Json::Value json;
  ifstream ifs(argv[1]);
  if ( ifs.fail() ) {
    cerr << "[ERROR] can not open " << argv[1] << endl;
    return __LINE__;
  }
  bool flag = reader.parse(ifs, json);
  if ( !flag ) {
    cerr << "[ERROR] " << reader.getFormattedErrorMessages() << endl;
    return __LINE__;
  }

  boost::filesystem::create_directories(json["outdir"].asString());
  char outfile[256];

  eigen_mati_t tris; eigen_matd_t nods;
  igl::readOBJ(json["mesh"].asString(), nods, tris);
  sprintf(outfile, "%s/origin.obj", json["outdir"].asString().c_str());
  igl::writeOBJ(outfile, nods, tris);

  boundary_loop loop;
  select_method(loop, static_cast<boundary_loop::type>(json["method"].asInt()));

  get_boundary_loop(tris, nods, loop);
  calc_bnd_local_frame(tris, nods, loop);
  calc_infinitesimal_elem(loop);
  calc_bnd_indicator(loop);
  cout << loop.sf << endl;

  eigen_vecd_t center(3);
  center << json["center"][0].asDouble(), json["center"][1].asDouble(), json["center"][2].asDouble();
  double dx = json["dx"].asDouble();
  eigen_matd_t bbox(2, 3), pts;
  bbox << center(0)-dx, center(1)-dx, center(2)-dx,
          center(0)+dx, center(1)+dx, center(2)+dx;
  cout << "[INFO] sampling bounding box:\n" << bbox << endl;
  structured_grid_sampling(bbox, json["resolution"].asInt(), pts);

  eigen_veci_t p = eigen_veci_t::Zero(pts.rows());
  int count = 0;
  std::generate(p.data(), p.data()+p.size(), [&count]()->int {return count++;});
  sprintf(outfile, "%s/grid.m%d.vtk", json["outdir"].asString().c_str(), json["method"].asInt());
  ofstream os(outfile);
  point2vtk(os, pts.data(), pts.rows(), p.data(), p.size());

  eigen_vecd_t sf;
  calc_scalar_field(loop, pts, sf);
  point_data(os, sf.data(), sf.size(), "sf");

  loop.T *= 0.1; loop.B *= 0.1; loop.N *= 0.1;
  string outx = json["outdir"].asString()+string("/T.vtk");
  draw_vert_direct_field(outx.c_str(), loop.pos.data(), loop.pos.rows(), loop.T.data());
  string outy = json["outdir"].asString()+string("/B.vtk");
  draw_vert_direct_field(outy.c_str(), loop.pos.data(), loop.pos.rows(), loop.B.data());
  string outz = json["outdir"].asString()+string("/N.vtk");
  draw_vert_direct_field(outz.c_str(), loop.pos.data(), loop.pos.rows(), loop.N.data());

  cout << "[INFO] done\n";
  return 0;
}
