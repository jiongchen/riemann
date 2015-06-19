#ifndef DEFORMATION_TRANSFER_H
#define DEFORMATION_TRANSFER_H

#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>

#include "def.h"

namespace geom_deform {

class deform_transfer
{
public:
  typedef zjucad::matrix::matrix<size_t> mati_t;
  typedef zjucad::matrix::matrix<double> matd_t;
  enum EnergyType {
    SMOOTHNESS,
    IDENTITY,
    CLOSEST
  };
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
  int solve_corres_precompute();
  int solve_corres_first_phase();
  int solve_corres_second_phase();
  // deformation solver
  int transfer_deformation();
  // debug
  int see_ghost_tet_mesh(const char *filename, const std::string &which) const;
  int see_corres_mesh(const char *filename) const;
private:
  void append_fourth_vert(const mati_t &tri_cell, const matd_t &tri_nods, mati_t &tet_cell, matd_t &tet_nods) const;
  void remove_fourth_vert(const mati_t &tet_cell, const matd_t &tet_nods, mati_t &tri_cell, matd_t &tri_nods) const;

  mati_t src_tris_, tar_tris_;
  matd_t src_ref_nods_, tar_ref_nods_;
  matd_t src_cor_nods_;
  matd_t src_def_nods_, tar_def_nods_;

  std::vector<std::tuple<size_t, size_t>> vert_map_;
  Eigen::MatrixXd Vinv_;

  std::vector<std::shared_ptr<surfparam::Functional<double>>> buff_;
  std::shared_ptr<surfparam::Functional<double>> energy_;
};

}
#endif
