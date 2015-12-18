#include "geometry_extend.h"

#include <Eigen/Eigen>
#include <jtflib/mesh/util.h>

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;
using mati_t=zjucad::matrix::matrix<size_t>;
using matd_t=zjucad::matrix::matrix<double>;

namespace riemann {

void gen_rand_orth_vec_3d(const double *x, double *xT) {
  Map<const Vector3d> X(x);
  Map<Vector3d> XT(xT);
  Vector3d rand_vec;
  do {
    rand_vec.setRandom();
    XT = X.cross(rand_vec);
  } while ( XT.norm() < 1e-8 );
}

int calc_vert_local_frame(const mati_t &tris, const matd_t &nods, matd_t &frame) {
  frame.resize(9, nods.size(2));
  {
    matd_t normal;
    jtf::mesh::cal_point_normal(tris, nods, normal);
    frame(colon(6, 8), colon()) = normal;
  }
#pragma omp parallel for
  for (size_t i = 0; i < nods.size(2); ++i) {
    frame(colon(6, 8), i) /= norm(frame(colon(6, 8), i));
    gen_rand_orth_vec_3d(&frame(6, i), &frame(0, i));
    frame(colon(0, 2), i) /= norm(frame(colon(0, 2), i));
    frame(colon(3, 5), i) = cross(frame(colon(6, 8), i), frame(colon(0, 2), i));
  }
  return 0;
}

void project_vector_on_subspace(const double *vect, const size_t dim,
                                const double *basis, const size_t sub_dim,
                                double *proj_vect) {
  Map<const VectorXd> b(vect, dim);
  Map<const MatrixXd> A(basis, dim, sub_dim);
  Map<VectorXd> pv(proj_vect, dim);
  MatrixXd LHS = A.transpose()*A;
  VectorXd rhs = A.transpose()*b;
  FullPivLU<MatrixXd> sol;
  sol.compute(LHS);
  VectorXd x = sol.solve(rhs);
  pv = A*x;
}

void project_point_on_subspace(const double *x, const size_t dim,
                               const double *origin, const double *basis, const size_t sub_dim,
                               double *proj_x) {
  Map<const VectorXd> X(x, dim), O(origin, dim);
  Map<VectorXd> Px(proj_x, dim);
  VectorXd OX = X-O, dx(dim);
  project_vector_on_subspace(&OX[0], dim, basis, sub_dim, &dx[0]);
  Px = O+dx;
}

void project_point_on_plane(const double *x, const double *origin,
                            const double *tan0, const double *tan1,
                            double *proj_x) {
  Map<const Vector3d> X(x), O(origin), T0(tan0), T1(tan1);
  Map<Vector3d> Px(proj_x);
  Matrix<double, 3, 2> A;
  A.col(0) = T0;
  A.col(1) = T1;
  Matrix2d ATA = A.transpose()*A;
  Vector2d b = A.transpose()*(X-O);
  Vector2d u(2);
  const double detATA = ATA(0,1)*ATA(1,0)-ATA(0,0)*ATA(1,1);
  u[0] = -(ATA(1,1)*b[0]-ATA(0,1)*b[1])/detATA;
  u[1] = (ATA(1,0)*b[0]-ATA(0,0)*b[1])/detATA;
  Px = O+A*u;
}

int calc_one_ring_face(const matrix<size_t> &tris, vector<vector<size_t>> &p2f) {
  if ( tris.size(1) != 3 ) {
    cerr << "[info] calc_one_ring_face: not triangle mesh!\n";
    return __LINE__;
  }
  const size_t vert_num = zjucad::matrix::max(tris)+1;
  p2f.resize(vert_num);
  for (size_t i = 0; i < tris.size(2); ++i) {
    p2f[tris(0, i)].push_back(i);
    p2f[tris(1, i)].push_back(i);
    p2f[tris(2, i)].push_back(i);
  }
  return 0;
}

}
