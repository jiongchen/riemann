#ifndef VOLUME_FRAME_H
#define VOLUME_FRAME_H

#include <zjucad/matrix/matrix.h>

#include "def.h"

namespace riemann {

static const double g_RXYZ [24][3][3]={
  {{1,0,0},{0,0,1},{0,-1,0}},  //  u1
  {{1,0,0},{0,-1,0},{0,0,-1}}, //  u2
  {{1,0,0},{0,0,-1},{0,1,0}},  //  u3
  {{0,0,-1},{0,1,0},{1,0,0}},  //  v1
  {{-1,0,0},{0,1,0},{0,0,-1}}, //  v2
  {{0,0,1},{0,1,0},{-1,0,0}},  //  v3
  {{0,1,0},{-1,0,0},{0,0,1}},  //  w1
  {{-1,0,0},{0,-1,0},{0,0,1}}, //  w2
  {{0,-1,0},{1,0,0},{0,0,1}},  //  w3
  {{1,0,0},{0,1,0},{0,0,1}},   //  I
  {{-1,0,0},{0,0,1},{0,1,0}},   //  u1 * v2, u3 * w2, w2 * u1
  {{-1,0,0},{0,0,-1},{0,-1,0}}, //  v2 * u1, u3 * v2
  {{0,0,1},{1,0,0},{0,1,0}},    //  v3 * u3, w3 * v3,  u3*w3
  {{0,0,-1},{1,0,0},{0,-1,0}},  //  v1 * u1, w3 * v1,
  {{0,0,-1},{-1,0,0},{0,1,0}},  //  v1 * u3, u3 * w1,
  {{0,0,1},{-1,0,0},{0,-1,0}},  //  v3 * u1, u1 * w1
  {{0,1,0},{0,0,1},{1,0,0}},    //  u1 * v1, v1 * w1, w1 * u1
  {{0,-1,0},{0,0,-1},{1,0,0}},  //  u3 * v1, w3 * u3
  {{0,0,1},{0,-1,0},{1,0,0}},   //  u2 * v1, v3 * u2
  {{0,1,0},{0,0,-1},{-1,0,0}},  //  u3 * v3, v3 * w1
  {{0,-1,0},{0,0,1},{-1,0,0}},  //  u1 * v3, w3 * u1, v3 * w3,
  {{0,0,-1},{0,-1,0},{-1,0,0}}, //  v1 * u2, w2 * v1, v3 * w2
  {{0,1,0},{1,0,0},{0,0,-1}},   //  u2 * w3, v2 * w1, w3 * v2
  {{0,-1,0},{-1,0,0},{0,0,-1}}  //  u2 * w1, w3 * u2
};

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

struct cross_frame_args
{
  double ws, wa;
  double epsf;
  size_t maxits;
};

class cross_frame_opt
{
public:
  static cross_frame_opt* create(const mati_t &tets, const matd_t &nods, const cross_frame_args &args);
  int solve_laplacian(Eigen::VectorXd &Fs) const;
  int solve_initial_frames(const Eigen::VectorXd &Fs, Eigen::VectorXd &abc) const;
  int optimize_frames(Eigen::VectorXd &abc) const;
private:
  int init(const mati_t &tets, const matd_t &nods, const cross_frame_args &args);
private:
  cross_frame_args args_;
  std::vector<std::shared_ptr<Functional<double>>> buffer_;
  std::shared_ptr<Functional<double>> energy_;
};

}

#endif
