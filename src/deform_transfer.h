#ifndef DEFORMATION_TRANSFER_H
#define DEFORMATION_TRANSFER_H

#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>

namespace geom_deform {

class deform_transfer
{
public:
  typedef zjucad::matrix::matrix<size_t> mati_t;
  typedef zjucad::matrix::matrix<double> matd_t;
  deform_transfer();
  // io
  int load_reference_source_mesh(const char *filename);
  int load_reference_target_mesh(const char *filename);
  int load_deformed_source_mesh(const char *filename);
  int load_vertex_markers(const char *filename);
  int save_reference_source_mesh(const char *filename) const;
  int save_reference_target_mesh(const char *filename) const;
  int save_deformed_source_mesh(const char *filename) const;
  int save_deformed_target_mesh(const char *filename) const;
  // correspondence solver
  // deformation solver
private:
  mati_t src_tris_, tar_tris_;
  matd_t src_ref_nods_, src_def_nods_;
  matd_t tar_ref_nods_, tar_def_nods_;
  Eigen::MatrixXd Vinv_;
};

}
#endif
