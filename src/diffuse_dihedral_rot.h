#ifndef DIFFUSE_DIHEDRAL_ROT_H
#define DIFFUSE_DIHEDRAL_ROT_H

#include <zjucad/matrix/matrix.h>

#include "def.h"

namespace riemann {

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

class graph_t;
typedef graph_t tree_t;

void diffuse_rotation(const mati_t &tris, const matd_t &restv, const matd_t &currv, const size_t root_face, const tree_t &g,
                      std::vector<Eigen::Matrix3d> &rot);

//class diffuse_arap : public Functional<double>
//{
//public:
//  diffuse_arap(const mati_t &tets, const matd_t &nods, const std::vector<Eigen::Matrix3d> &R);
//  size_t Nx() const;
//  int Val(const double *x, double *val) const;
//  int Gra(const double *x, double *gra) const;
//  int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
//  void pin_down_verts(const size_t id, const double *pos);
//private:
//  std::vector<Eigen::Matrix3d> R_; //face size
//};

}

#endif
