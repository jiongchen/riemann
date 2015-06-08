#ifndef FRAME_FIELD_DEFORM_H
#define FRAME_FIELD_DEFORM_H

#include <Eigen/Sparse>
#include <zjucad/matrix/matrix.h>

namespace geom_deform {

class frame_field_deform
{
public:
  typedef zjucad::matrix::matrix<size_t> mati_t;
  typedef zjucad::matrix::matrix<double> matd_t;
  frame_field_deform();
  // IO
  int load_mesh(const char *file);
  int load_constraints(const char *file);
  int save_original_mesh(const char *file) const;
  int save_deformed_mesh(const char *file) const;
  // prepare
  int interp_frame_fields();
  // deform
  int precompute();
  int deform();
  // debug
  int save_local_frame(const char *file, const double len) const;
private:
  int build_local_frames();
  mati_t tris_;
  matd_t nods_, _nods_;
  std::vector<Eigen::Matrix3d> B_;
  std::vector<size_t> fidx_;
  std::vector<size_t> g2l_;
  Eigen::VectorXd W_;
};

}
#endif
