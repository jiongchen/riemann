#ifndef FRAME_FIELD_DEFORM_H
#define FRAME_FIELD_DEFORM_H

#include <Eigen/Sparse>

namespace geom_deform {

class frame_field_deform
{
public:
  typedef Eigen::MatrixXd matd_t;
  typedef Eigen::MatrixXi mati_t;
  typedef Eigen::VectorXd vec_t;
  typedef Eigen::SparseMatrix<double> spmat_t;
  frame_field_deform();
  // IO
  int load_mesh(const char *path);
  int load_constraints(const char *path);
  int save_origin_mesh(const char *path) const;
  int save_deformed_mesh(const char *path) const;

  // prepare
  int interp_frame_fields();

  // deform
  int deform();

private:
  mati_t tris_;
  matd_t nods_;
  std::vector<Eigen::Matrix2d> B_;
  std::vector<Eigen::Vector3d> W_;
};

}
#endif
