#ifndef GRADIENT_BASED_DEFORM_H
#define GRADIENT_BASED_DEFORM_H

#include <Eigen/Sparse>
#include <zjucad/matrix/matrix.h>
#include <unordered_set>

namespace riemann {

/// Lp = div w;
/// distinguish the gradient to coordinates
/// on surface and the deformation gradient
class gradient_field_deform
{
public:
  using mati_t = zjucad::matrix::matrix<size_t>;
  using matd_t = zjucad::matrix::matrix<double>;
  gradient_field_deform();
  // io
  int load_origin_model(const char *filename);
  int save_origin_model(const char *filename) const;
  int save_deformed_model(const char *filename) const;
  // init
  int init();
  // manipulate
  int scale_grad_fields(const double scale);
  int rotate_grad_fields(const double *axis, const double angle);
  int reverse_grad_fields();
  int set_fixed_verts(const std::vector<size_t> &idx);
  int manipualte(const size_t idx, const double *u);
  // precompute
  int precompute();
  // deform
  int deform();
  // debug
  int see_coord_grad_fields(const char *filename, const int xyz) const;
  int see_harmonic_field(const char *filename) const;
private:
  int calc_init_coord_grad();
  int calc_bary_basis_grad();
  int calc_element_area();
  int solve_for_xyz(const int xyz);
  int calc_divergence(const Eigen::VectorXd &vf, Eigen::VectorXd &div) const;

  mati_t tris_;
  matd_t nods_, _nods_;
  Eigen::SparseMatrix<double> G_;   // 3*#face by #vert
  Eigen::SparseMatrix<double> L_;   // #vert by #vert
  Eigen::MatrixXd grad_xyz_;
  Eigen::MatrixXd gradB_;           // 3 by 3*#face
  Eigen::VectorXd area_;

  std::unordered_set<size_t> fixDOF_;
  std::vector<size_t> g2l_;
  Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> sol_;
  Eigen::VectorXd hf_;
};

}
#endif
