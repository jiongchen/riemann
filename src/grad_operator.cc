#include "grad_operator.h"

#include <zjucad/matrix/itr_matrix.h>

#include "geometry_extend.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

namespace riemann {

void calc_grad_operator(const mati_t &tris, const matd_t &nods, spmat_t *G) {
  vector<Triplet<double>> trips;
  /// for every face
  for (size_t i = 0; i < tris.size(2); ++i) {
    matd_t vert = nods(colon(), tris(colon(), i));
    matd_t H(3, 3);
    calc_tri_height_vector<3>(&vert[0], &H[0]);
    for (size_t j = 0; j < 3; ++j) {
      double w = dot(H(colon(), j), H(colon(), j));
      H(colon(), j) /= w;
    }
    for (size_t p = 0; p < 3; ++p)
      for (size_t q = 0; q < 3; ++q)
        trips.push_back(Triplet<double>(3*i+p, tris(q, i), H(p, q)));
  }
  G->resize(3*tris.size(2), nods.size(2));
  G->reserve(trips.size());
  G->setFromTriplets(trips.begin(), trips.end());
}

void calc_tet_height_vector(const double *vert, double *height) {
  itr_matrix<const double *> X(3, 4, vert);
  itr_matrix<double *> H(3, 4, height);
  for (size_t i = 0; i < 4; ++i)  {
    matd_t projX(3, 1);
    matd_t u = X(colon(), (i+2)%4)-X(colon(), (i+1)%4);
    matd_t v = X(colon(), (i+3)%4)-X(colon(), (i+1)%4);
    project_point_on_plane(&X(0, i), &X(0, (i+1)%4), &u[0], &v[0], &projX[0]);
    H(colon(), i) = X(colon(), i)-projX;
  }
}

void calc_tet_linear_basis_grad(const double *vert, double *gradB) {
  calc_tet_height_vector(vert, gradB);
  Map<Matrix<double, 3, 4>> G(gradB);
  G.col(0) /= G.col(0).squaredNorm();
  G.col(1) /= G.col(1).squaredNorm();
  G.col(2) /= G.col(2).squaredNorm();
  G.col(3) /= G.col(3).squaredNorm();
}

}
