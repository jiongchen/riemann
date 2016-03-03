#include "hole_filling.h"

#include <random>
#include <chrono>
#include <Eigen/UmfPackSupport>

#include "igl/boundary_loop.h"
#include "igl/per_vertex_normals.h"
#include "igl/cross.h"
#include "config.h"

using namespace std;
using namespace Eigen;

namespace riemann {

void get_boundary_loop(const eigen_mati_t &tris, const eigen_matd_t &nods, boundary_loop &loop) {
  igl::boundary_loop(tris, loop.bnd);
  loop.pos.resize(loop.bnd.size(), 3);
  for (int i = 0; i < loop.pos.rows(); ++i)
    loop.pos.row(i) = nods.row(loop.bnd(i));
}

void select_method(boundary_loop &loop, boundary_loop::type t) {
  loop.method = static_cast<int>(t);
  if ( loop.method != 0 && loop.method != 1 ) {
    cerr << "[ERROR] method not supported!\n";
    exit(EXIT_FAILURE);
  }
}

void calc_bnd_local_frame(const eigen_mati_t &tris, const eigen_matd_t &nods, boundary_loop &loop) {
  eigen_matd_t vn;
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

void calc_infinitesimal_elem(boundary_loop &loop) {
  const int bnd_size = loop.bnd.size();
  loop.dl.setZero(bnd_size);
  for (int i = 0; i < bnd_size; ++i) {
    loop.dl(i) += (loop.pos.row(i)-loop.pos.row((i-1+bnd_size)%bnd_size)).norm();
    loop.dl(i) += (loop.pos.row(i)-loop.pos.row((i+1)%bnd_size)).norm();
    loop.dl(i) /= 2.0;
  }
}

void calc_bnd_indicator(boundary_loop &loop) {
  const int bnd_size = loop.bnd.size();
  loop.sf.setZero(bnd_size);
  if ( loop.method == 0 ) {
    for (int i = 0; i < bnd_size; ++i) {
      double integral = 0.0;
      for (int k = 0; k < bnd_size; ++k) {
        integral += (loop.pos.row(i)-loop.pos.row(k)).dot(loop.B.row(k))*loop.dl(k);
      }
      loop.sf(i) = -integral;
    }
  } else {
    MatrixXd A(bnd_size, bnd_size); VectorXd b(bnd_size);
    b.setZero();
    for (int i = 0; i < bnd_size; ++i) {
      for (int j = 0; j < bnd_size; ++j) {
        A(i, j) = (loop.pos.row(i)-loop.pos.row(j)).dot(loop.B.row(j))*loop.dl(j);
      }
    }
    cout << A.rows() << endl;
    FullPivLU<MatrixXd> solver;
    solver.compute(A);
    cout << solver.rank() << endl;
    getchar();
    VectorXd x = solver.solve(b);
    std::copy(x.data(), x.data()+x.size(), loop.sf.data());
  }
}

static inline double calc_rbf(const eigen_vecd_t &x, const eigen_vecd_t &c) {
  return 1.0/(x-c).norm();
}

double indicator_value(const eigen_vecd_t &x, const boundary_loop &loop) {
  double rtn = 0.0;
  if ( loop.method == 0 ) {
    double integral = 0.0, correction = 0.0, sum = 0.0;
    // first part
    for (int i = 0; i < loop.bnd.size(); ++i) {
      integral += ((x-loop.pos.row(i)).dot(loop.B.row(i)))*loop.dl(i);
    }
    // second part
    for (int i = 0; i < loop.bnd.size(); ++i) {
      double rbf = calc_rbf(x, loop.pos.row(i));
      correction += rbf*loop.sf(i)*loop.dl(i);
      sum += rbf*loop.dl(i);
    }
    rtn = integral+correction/sum;
  } else {
    for (int i = 0; i < loop.bnd.size(); ++i)
      rtn += ((x-loop.pos.row(i)).dot(loop.sf(i)*loop.B.row(i)))*loop.dl(i);
  }
  return rtn;
}

void uniform_random_sampling(const eigen_matd_t &box, const int num, eigen_matd_t &pts) {
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

void structured_grid_sampling(const eigen_matd_t &box, const int res, eigen_matd_t &pts) {
  double xmin = box(0, 0), xmax = box(1, 0);
  double ymin = box(0, 1), ymax = box(1, 1);
  double zmin = box(0, 2), zmax = box(1, 2);
  vector<double> buffer;
  double dx = (xmax-xmin)/res, dy = (ymax-ymin)/res, dz = (zmax-zmin)/res;
  for (int i = 0; i <= res; ++i) {
    for (int j = 0; j <= res; ++j) {
      for (int k = 0; k <= res; ++k) {
        buffer.push_back(xmin+i*dx);
        buffer.push_back(ymin+j*dy);
        buffer.push_back(zmin+k*dz);
      }
    }
  }
  pts.resize(buffer.size()/3, 3);
  std::copy(buffer.begin(), buffer.end(), pts.data());
}

void calc_scalar_field(const boundary_loop &loop, const eigen_matd_t &pts, eigen_vecd_t &sf) {
  sf.setZero(pts.rows());
#pragma omp parallel for
  for (int i = 0; i < sf.size(); ++i) {
    sf(i) = indicator_value(pts.row(i), loop);
  }
}

}
