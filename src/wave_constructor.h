#ifndef WAVE_CONSTRUCTOR_H
#define WAVE_CONSTRUCTOR_H

#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>
#include <Eigen/UmfPackSupport>
#include <jtflib/mesh/mesh.h>

namespace riemann {

template <typename T>
class Constraint;

class wave_constructor
{
public:
  using mati_t = zjucad::matrix::matrix<size_t>;
  using matd_t = zjucad::matrix::matrix<double>;
  // io
  int load_model_from_obj(const char *filename);
  int load_frame_field(const char *filename);
  int load_feature_line(const char *filename);
  int save_feature_to_vtk(const char *filename) const;
  int save_model_to_obj(const char *filename) const;
  int save_wave_to_vtk(const char *filename) const;
  // solve
  int scale_frame_field(const double scale);
  int solve_phase_transition();
  int prepare();
  int give_an_initial_value(const size_t idx);
  int solve_wave();
  // debug
  int vis_edge_frame_field(const char *file_x, const char *file_y, const double scale) const;
private:
  int extract_edges();
  int build_frame_on_vert();
  void count_vert_show_on_feature();
  int solve_wave_soft_feature();
  int solve_wave_hard_feature();

  mati_t tris_;
  matd_t nods_;
  mati_t lines_;
  matd_t vert_frm_;
  matd_t edge_frm_;  // 6 x #edge
  mati_t edges_;     // 2 x #edge
  std::map<size_t, size_t> vert_count_;

  matd_t f_;         // 4 x #vert
  matd_t cIJ_, cJI_; // 4 x #edge

  std::unique_ptr<jtf::mesh::edge2cell_adjacent> e2c_;
  std::vector<std::shared_ptr<Constraint<double>>> buff_;
  std::shared_ptr<Constraint<double>> constraint_;
  std::shared_ptr<Constraint<double>> feature_cons_;
  Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> ltl_solver_;
  Eigen::UmfPackLU<Eigen::SparseMatrix<double>> lu_solver_;
};

}

#endif
