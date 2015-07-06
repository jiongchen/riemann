#include "grad_operator.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

namespace geom_deform {

/// 0------->1
///          |
///          |
///          v
///          2
void calc_tri_height_vector(const double *vert, double *height) {
  Map<const Matrix3d> X(vert);
  Map<Matrix3d> H(height);
  Matrix3d I = Matrix3d::Identity();
  for (size_t i = 0; i < 3; ++i) {
    size_t idx[3] = {i, (i+1)%3, (i+2)%3};
    Vector3d x12 = X.col(idx[2])-X.col(idx[1]);
    H.col(i) = (I-x12*x12.transpose()/x12.dot(x12))*(X.col(idx[0])-X.col(idx[1]));
  }
}

void calc_grad_operator(const mati_t &tris, const matd_t &nods, spmat_t *G) {
  vector<Triplet<double>> trips;
  /// for every face
  for (size_t i = 0; i < tris.size(2); ++i) {
    matd_t vert = nods(colon(), tris(colon(), i));
    matd_t H(3, 3);
    calc_tri_height_vector(&vert[0], &H[0]);
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

}
