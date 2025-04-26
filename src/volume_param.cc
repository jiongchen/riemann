#include "volume_param.h"

#include "def.h"
#include "blas_lapack.h"
#include <zjucad/matrix/lapack.h>

namespace riemann {

// class volume_param_energy : public Functional<double>
// {
//  public:
//   volume_param_energy(const mati_t &tets, const matd_t &nods)
//       : tets_(tets), dim_(nods.size()) {
//     vol_.resize(tets.size(2));
//     Dm_.resize(9, tets.size(2));
//     #pragma omp parallel for
//     for (size_t i = 0; i < tets_.size(2); ++i) {
//       matd_t Ds = nods(colon(), tets(colon(1, 3), i))-nods(colon(), tets(0, i))*ones<double>(1, 3);
//       matd_t Ds_copy = Ds;
//       vol_[i] = fabs(det(Ds_copy))/6.0;
//       inv(Ds);
//       Dm_(colon(), i) = Ds(colon());
//     }
//     RT_.resize(9, tets.size(2));
//   }
//   size_t Nx() const {
//     return dim_;
//   }
//   int Val(const double *x, double *val) const {
//     itr_matrix<const double *> X(3, dim_/3, x);
//     matd_t vert = zeros<double>(3, 4); double value = 0;
//     for (size_t i = 0; i < tets_.size(2); ++i) {
//       vert = X(colon(), tets_(colon(), i));
//     }
//     return 0;
//   }
//   int Gra(const double *x, double *gra) const {
//     return 0;
//   }
//   int Hes(const double *x, vector<Triplet<double>> *hes) const {
//     return 0;
//   }
//   void set_frames();
//  private:
//   const mati_t &tets_;
//   matd_t Dm_;
//   matd_t vol_;
//   matd_t RT_;
// };

// volume_param::volume_param(const mati_t &tets, const matd_t &nods)
//     : tets_(tets) {
//   arap_ = make_shared<>();
// }

// int volume_param::comb_frames(double *frames) const {
//   return 0;
// }

// int volume_param::parameterize(const double *frames, double *x) const {
//   return 0;
// }

}
