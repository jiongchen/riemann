#ifndef WAVE_CONSTRUCTOR_H
#define WAVE_CONSTRUCTOR_H

#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>
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
  int save_model_to_obj(const char *filename) const;
  int save_wave_to_vtk(const char *filename) const;
  // solve
  int init();
  int solve_phase_transition();
  int prepare();
  int solve_wave_value();
  // debug
  void test_wave_conditions() const;
  int vis_edge_frame_field(const char *file_x, const char *file_y, const double scale) const;
private:
  int extract_edges();
  int build_frame_on_vert();
  void count_vert_on_features();

  mati_t tris_;
  matd_t nods_;
  matd_t local_frame_;
  matd_t frame_field_;   // 6 x #edge
  matd_t size_field_;    // 2 x #edge
  mati_t edges_;         // 2 x #edge
  std::map<size_t, size_t> vert_on_fl_;

  matd_t f_;         // 4 x #vert
  matd_t alphaIJ_;   // 1 x #edge
  matd_t betaIJ_;    // 1 x #edge
  matd_t cIJ_, cJI_; // 4 x #edge

  std::unique_ptr<jtf::mesh::edge2cell_adjacent> e2c_;
  std::vector<std::shared_ptr<Constraint<double>>> buff_;
  std::shared_ptr<Constraint<double>> constraint_;
  Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver_;
};

}

#endif
