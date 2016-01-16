#ifndef GRADIENT_BASED_DEFORM_H
#define GRADIENT_BASED_DEFORM_H

#include <Eigen/Sparse>
#include <zjucad/matrix/matrix.h>
#include <unordered_set>

namespace riemann {

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

/// Lp = div w;
/// distinguish the gradient to coordinates
/// on surface and the deformation gradient
class gradient_field_deform
{
public:
  gradient_field_deform(const mati_t &tris, const matd_t &nods);
  void set_fixed_verts(const std::vector<size_t> &idx);
  void set_edited_verts(const std::vector<size_t> &idx);
//  int edit_boundary(const std::vector<size_t> &idx, uniform transform);
  int solve_harmonic_field(); /// todo
  int deform(double *x);
public:
  Eigen::MatrixXd Gxyz_;
  Eigen::VectorXd hf_;
private:
  int calc_bary_basis_grad(const mati_t &tris, const matd_t &nods, Eigen::MatrixXd &gradB);
  int solve_xyz(const int xyz);
  int calc_div(const Eigen::VectorXd &vf, Eigen::VectorXd &div) const;
  int precompute();
private:
  const mati_t &tris_;
  const matd_t &nods_;

  Eigen::SparseMatrix<double> G_, L_;
  Eigen::MatrixXd gradB_;           // 3 by 3*#face
  Eigen::VectorXd area_;

  Eigen::MatrixXd transform_;       // 5 by #vert
  std::unordered_set<size_t> fixDoF_, editDoF_;
  std::vector<size_t> g2l_;
};

}
#endif
