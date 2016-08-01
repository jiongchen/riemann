#include "volume_param.h"

#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>

#include "def.h"

namespace riemann {

class volume_param_energy : public Functional<double>
{
};

volume_param::volume_param(const mati_t &tets, const matd_t &nods)
    : tets_(tets) {
  vol_.resize(tets.size(2));
  Dm_.resize(9, tets.size(2));
  #pragma omp parallel for
  for (size_t i = 0; i < tets_.size(2); ++i) {
    matd_t Ds = nods(colon(), tets(colon(1, 3), i))-nods(colon(), tets(0, i))*ones<double>(1, 3);
    matd_t Ds_copy = Ds;
    vol_[i] = fabs(det(Ds_copy))/6.0;
    inv(Ds);
    Dm_(colon(), i) = Ds(colon());
  }
}

int volume_param::comb_frames(double *frames) const {
  return 0;
}

int volume_param::parameterize(const double *frames, double *x) const {
  return 0;
}

}
