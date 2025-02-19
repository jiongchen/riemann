#ifndef VOLUME_FRAME_H
#define VOLUME_FRAME_H

#include <Eigen/Dense>
#include <zjucad/matrix/matrix.h>
#include <boost/property_tree/ptree.hpp>

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
using boost::property_tree::ptree;

template <typename T>
class Functional;

void convert_zyz_to_mat(const Eigen::VectorXd &abc, Eigen::VectorXd &mat);

class cross_frame_opt
{
public:
  cross_frame_opt(const mati_t &tets, const matd_t &nods, const ptree &pt);
  int solve_laplacian(Eigen::VectorXd &Fs) const;
  int solve_initial_frames(const Eigen::VectorXd &Fs, Eigen::VectorXd &abc) const;
  int optimize_frames(Eigen::VectorXd &abc) const;
private:
  const mati_t &tets_;
  const matd_t &nods_;
  const ptree &pt_;
  std::vector<std::shared_ptr<Functional<double>>> buffer_;
};

class frame_smoother
{
public:
  frame_smoother(const mati_t &tets, const matd_t &nods, const ptree &pt);
  int smoothSH(Eigen::VectorXd &abc) const;
  int smoothL1(Eigen::VectorXd &mat) const;
private:
  const mati_t &tets_;
  const matd_t &nods_;
  const ptree &pt_;
  zjucad::matrix::matrix<int> is_bnd_tet_;
};

}

#endif
