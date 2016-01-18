#ifndef GRADIENT_BASED_DEFORM_H
#define GRADIENT_BASED_DEFORM_H

#include <Eigen/Sparse>
#include <zjucad/matrix/matrix.h>
#include <unordered_set>
#include <Eigen/Geometry>

namespace riemann {

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

/// Lp = div w;
class gradient_field_deform
{
public:
  gradient_field_deform(const mati_t &tris, const matd_t &nods);
  void set_fixed_verts(const std::vector<size_t> &idx);
  void set_edited_verts(const std::vector<size_t> &idx);
  void prescribe_uniform_transform(const Eigen::Quaterniond &q, const double s); // local rigid
  int precompute();
  int deform(double *x) const;
  // internal in future
    int solve_harmonic_field();
    void propogate_transform();
    void interp_face_transform();
    void rotate_surf_piecewise(mati_t &rtris, matd_t &rnods) const;
    //void update_gradient_field(const mati_t &rtris, const matd_t &rnods);
    void update_gradient_field();
    void prefactorize();
private:
  int solve_xyz(double *x, const int xyz) const;
  int calc_div(const Eigen::VectorXd &vf, Eigen::VectorXd &div) const;
public:
  const size_t dim_;
  const mati_t &tris_;
  const matd_t &nods_;
  std::unordered_set<size_t> fixDoF_, editDoF_;

  Eigen::SparseMatrix<double> L_;
  Eigen::MatrixXd gradShape_;
  Eigen::VectorXd area_;
  Eigen::VectorXd hf_;

  const Eigen::Quaterniond ID_;
  Eigen::Quaterniond q_;
  double s_;
  std::vector<Eigen::Quaterniond> vR_, fR_;
  Eigen::VectorXd vs_, fs_;

  Eigen::MatrixXd Gxyz_;
  std::vector<size_t> g2l_;
  Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver_;
};

}
#endif
