#ifndef WAVE_CONSTRUCTOR_H
#define WAVE_CONSTRUCTOR_H

#include <zjucad/matrix/matrix.h>

namespace riemann {

template <typename T>
class Functional;

class wave_constructor
{
public:
  using mati_t = zjucad::matrix::matrix<size_t>;
  using matd_t = zjucad::matrix::matrix<double>;
  // io
  int load_model_from_obj(const char *filename);
  int load_frame_field(const char *filename);
  int save_model_to_obj(const char *filename) const;
  int save_scalar_field_to_vtk(const char *filename) const;
  // solve
  int solve_phase_transition();
  int solve_wave_value();
  // debug
  void test_wave_conditions() const;
private:
  mati_t tris_;
  matd_t nods_;
  mati_t edge_;               // 2 x #edge
  matd_t frame_field_; // 3 x 2#edge
  matd_t size_field_;     // 1 x 2#edge

  matd_t f_; // 4 x #vert
  matd_t alpha_ij_; // 1 x #edge
  matd_t beta_ij_; // 1 x #edge
  matd_t c_ij_; // 4 x #edge

  std::vector<std::shared_ptr<Functional<double>>> buff_;
  std::shared_ptr<Functional<double>> energy_;
};

}

#endif
