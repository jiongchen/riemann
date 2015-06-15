#ifndef FRAME_FIELD_DEFORM_H
#define FRAME_FIELD_DEFORM_H

#include <Eigen/Sparse>
#include <zjucad/matrix/matrix.h>
#include <unordered_set>
#include <boost/property_tree/ptree.hpp>

#include "def.h"

namespace geom_deform {

class frame_field_deform
{
public:
  typedef zjucad::matrix::matrix<size_t> mati_t;
  typedef zjucad::matrix::matrix<double> matd_t;
  enum EnergyType {
    DEFORM,
    SMOOTH
  };
  frame_field_deform();
  frame_field_deform(const boost::property_tree::ptree &pt);
  // io
  int load_mesh(const char *file);
  int load_constraints(const char *file);
  int save_original_mesh(const char *file) const;
  int save_deformed_mesh(const char *file) const;
  int save_corss_field(const char *file) const;
  // prepare
  int interp_frame_fields();
  // deform
  int precompute();
  int deform();
  // frame to cross
  int gen_cross_field();
  // debug
  int visualize_local_bases(const char *file, const double len) const;
  int visualize_init_frames(const char *file, const double scale=1.0) const;
  int visualize_frame_fields(const char *file, const double scale=1.0);
  int visualize_tensor_fields(const char *file);
  bool check_spd_tensor_fields() const;
private:
  int build_local_bases(const mati_t &tris, const matd_t &nods, Eigen::MatrixXd &B);
  int interp_cross_fields();
  mati_t tris_;
  matd_t nods_, _nods_;
  Eigen::MatrixXd B_, _B_;

  std::vector<size_t> cons_face_;
  std::unordered_set<size_t> ffc_;
  std::vector<size_t> g2l_;
  Eigen::VectorXd W_;
  Eigen::MatrixXd X_;
  Eigen::MatrixXd F_;

  std::vector<std::shared_ptr<surfparam::Functional<double>>> buff_;
  std::shared_ptr<surfparam::Functional<double>> e_;
  Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> sol_;
  Eigen::SparseMatrix<double> LHS_;
  size_t max_iter_;
  double tolerance_;
  double lambda_;
  double perturb_;
};

}
#endif
