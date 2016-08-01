#ifndef VOLUME_PARAM_H
#define VOLUME_PARAM_H

#include <zjucad/matrix/matrix.h>

namespace riemann {

using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

template <typename T>
class Functional;

class volume_param
{
public:
  volume_param(const mati_t &tets, const matd_t &nods);
  int comb_frames(double *frames) const;
  int parameterize(const double *frames, double *x) const;
private:
  const mati_t &tets_;
  matd_t Dm_;
  matd_t vol_;
  std::shared_ptr<Functional<double>> arap_;
};

}

#endif
