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

class SH_smooth_energy : public Functional<double>
{
public:
  SH_smooth_energy(const mati_t &tets, const matd_t &nods, const double w);
  size_t Nx() const;
  int Val(const double *abc, double *val) const;
  int Gra(const double *abc, double *gra) const;
  int Hes(const double *abc, std::vector<Eigen::Triplet<double>> *hes) const { return __LINE__; }
private:
  double w_;
  size_t dim_;
  const mati_t &tets_;
  const matd_t &nods_;
  mati_t tet_face_;
  matd_t CR_;
  matd_t vol_;
};

class SH_align_energy : public Functional<double>
{
public:
  SH_align_energy(const mati_t &tets, const matd_t &nods, const double w);
  size_t Nx() const;
  int Val(const double *abc, double *val) const;
  int Gra(const double *abc, double *gra) const;
  int Hes(const double *abc, std::vector<Eigen::Triplet<double>> *hes) const { return __LINE__; }
private:
  double w_;
  size_t dim_;
  const mati_t &tets_;
  const matd_t &nods_;
  mati_t bnd_face_;
  matd_t zyz_;
  matd_t area_;
};

// class SH_project_energy
// {
// public:
//   SH_project_energy(const matd_t &f);
//   size_t Nx() const;
//   int Val(const double *abc, double *val) const;
//   int Gra(const double *abc, double *gra) const;
//   int Hes(const double *abc, double *hes) const;
// private:
  
// };

// class cross_frame_opt
// {
// public:
//   void solve_smooth_sh_coef();
//   void solve_initial_zyz();
//   void optimize_cross_frame();
// };

}

#endif
