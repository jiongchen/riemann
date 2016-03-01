#include "hole_filling.h"

#include <random>
#include <chrono>

#include "igl/boundary_loop.h"
#include "igl/per_vertex_normals.h"
#include "igl/cross.h"

using namespace std;
using namespace Eigen;

namespace riemann {

void get_boundary_loop(const mati &tris, const matd &nods, boundary_loop &loop) {
  igl::boundary_loop(tris, loop.bnd);
  loop.pos.resize(loop.bnd.size(), 3);
  for (int i = 0; i < loop.pos.rows(); ++i)
    loop.pos.row(i) = nods.row(loop.bnd(i));
}

void calc_bnd_local_frame(const mati &tris, const matd &nods, boundary_loop &loop) {
  matd vn;
  igl::per_vertex_normals(nods, tris, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA, vn);
  const int bnd_size = loop.bnd.size();
  loop.T.resize(bnd_size, 3);
  loop.B.resize(bnd_size, 3);
  loop.N.resize(bnd_size, 3);
  for (int i = 0; i < bnd_size; ++i) {
    loop.T.row(i) = (loop.pos.row((i+1)%bnd_size)-loop.pos.row(i)).normalized();
    loop.N.row(i) = vn.row(loop.bnd(i)).normalized();
    igl::cross(&loop.T(i, 0), &loop.N(i, 0), &loop.B(i, 0));
    loop.B.row(i) /= loop.B.row(i).norm();
  }
}

void calc_dl(boundary_loop &loop) {
  const int bnd_size = loop.bnd.size();
  loop.dl.setZero(bnd_size);
  for (int i = 0; i < bnd_size; ++i) {
    loop.dl(i) += (loop.pos.row(i)-loop.pos.row((i-1)%bnd_size)).norm();
    loop.dl(i) += (loop.pos.row(i)-loop.pos.row((i+1)%bnd_size)).norm();
    loop.dl(i) /= 2.0;
  }
}

//void init_bnd_scalar_field(boundary_loop &loop) {
//  const int bnd_size = loop.bnd.size();
//  loop.sf.setZero(bnd_size);
//  for (int i = 0; i < bnd_size; ++i) {
//    double integral = 0.0;
//    for (int k = 0; k < bnd_size; ++k) {
//      integral += (loop.pos.row(i)-loop.pos.row(k)).dot(loop.B.row(k))*loop.dl(k);
//    }
//    loop.sf(i) = integral/loop.dl.sum();
//  }
//}

//double indicator_value(const vecd &x, const boundary_loop &loop) {
//  double rtn = 0.0;
//  for (int i = 0; i < loop.bnd.size(); ++i) {
//    rtn += ((x-loop.pos.row(i)).dot(loop.B.row(i))-loop.sf(i))*loop.dl(i);
//  }
//  return rtn;
//}

void uniform_sampling(const matd &box, const int num, matd &pts) {
  double xmin = box(0, 0), xmax = box(1, 0);
  double ymin = box(0, 1), ymax = box(1, 1);
  double zmin = box(0, 2), zmax = box(1, 2);
  pts.resize(num, 3);
  unsigned seed = chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator(seed);
  uniform_real_distribution<double> distribution(0.0, 1.0);
  for (int i = 0; i < num; ++i) {
    pts(i, 0) = xmin+distribution(generator)*(xmax-xmin);
    pts(i, 1) = ymin+distribution(generator)*(ymax-ymin);
    pts(i, 2) = zmin+distribution(generator)*(zmax-zmin);
  }
}

}
