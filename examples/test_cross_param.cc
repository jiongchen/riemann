#include <iostream>
#include <fstream>
#include <zjucad/matrix/matrix.h>
#include <zjucad/matrix/itr_matrix.h>
#include <zjucad/matrix/io.h>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <boost/filesystem.hpp>

#include "src/vtk.h"
#include "src/grad_operator.h"
#include "src/util.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

void calc_barycentric_grad(const double *vert, double *height) {
  geom_deform::calc_tri_height_vector<2>(vert, height);
  Map<Matrix<double, 2, 3>> H(height);
  for (size_t j = 0; j < 3; ++j)
    H.col(j) /= H.col(j).squaredNorm();
}

class deform_energy
{
public:
  using mati_t = matrix<size_t>;
  using matd_t = matrix<double>;
  deform_energy(const mati_t &tris, const matd_t &nods, const Matrix2d &cross, const double h=1.0)
    : tris_(tris), nods_(nods), cross_(cross), h_(h) {
    calc_barycentric_grad(nods_.begin(), gradB_.data());
  }
  size_t dim() const {
    return nods_.size();
  }
  int val(const double *x, double *val) const {
    Map<const Matrix<double, 2, 3>> X(x);
    Vector3d u = X.row(0).transpose(), v = X.row(1).transpose();
    *val += (h_*gradB_*u-cross_.col(0)).squaredNorm()
        + (h_*gradB_*v-cross_.col(1)).squaredNorm();
    return 0;
  }
  int gra(const double *x, double *gra) const {
    Map<const Matrix<double, 2, 3>> X(x);
    Map<Matrix<double, 2, 3>> G(gra);
    Vector3d u = X.row(0).transpose(), v = X.row(1).transpose();
    G.row(0) += ((h_*gradB_).transpose()*(h_*gradB_*u-cross_.col(0))).transpose();
    G.row(1) += ((h_*gradB_).transpose()*(h_*gradB_*v-cross_.col(1))).transpose();
    return 0;
  }
  int hes(const double *x, double *hes) const {
    Matrix3d hess = h_*h_*gradB_.transpose()*gradB_;
    Map<Matrix<double, 6, 6>> H(hes);
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        H(2*i, 2*j) += hess(i, j);
        H(2*i+1, 2*j+1) += hess(i, j);
      }
    }
    return 0;
  }
private:
  const mati_t &tris_;
  const matd_t &nods_;
  const Matrix2d &cross_;
  const double h_;
  Matrix<double, 2, 3> gradB_;
};

#define QUERY_ENERGY_VALUE(e, x)                  \
  do {                                            \
    double v = 0;                                 \
    e->val(&x[0], &v);                            \
    cout << "[info] energy value: " << v << "\n"; \
  } while (0);

class cross_field_param
{
public:
  using mati_t = matrix<size_t>;
  using matd_t = matrix<double>;
  cross_field_param() {}
  // io
  int set_triangle() {
    srand(time(NULL));
    tris_ = colon(0, 2);
    nods_ = rand(2, 3);
    return 0;
  }
  int set_cross() {
    srand(time(NULL));
    cross_.col(0) = Vector2d::Random();
    cross_.col(0) /= cross_.col(0).norm();
    Matrix2d rot = Matrix2d::Zero();
    rot(0, 1) = -1;
    rot(1, 0) = 1;
    cross_.col(1) = rot*cross_.col(0);
    return 0;
  }
  int see_triangle(const char *filename) const {
    matd_t vert = zeros<double>(3, 3);
    vert(colon(0, 1), colon()) = nods_;
    ofstream os(filename);
    tri2vtk(os, &vert[0], vert.size(2), &tris_[0], tris_.size(2));
    os.close();
    return 0;
  }
  int see_cross(const char *filename) const {
    mati_t line = zeros<size_t>(2, 2);
    line(1, colon()) = colon(1, 2);
    matd_t vert(3, 3);
    vert(colon(0, 1), 0) = nods_*ones<double>(3, 1)/3.0;
    vert(colon(0, 1), 1) = vert(colon(0, 1), 0)+itr_matrix<const double*>(2, 1, &cross_(0, 0));
    vert(colon(0, 1), 2) = vert(colon(0, 1), 0)+itr_matrix<const double*>(2, 1, &cross_(0, 1));
    ofstream os(filename);
    line2vtk(os, &vert[0], vert.size(2), &line[0], line.size(2));
    os.close();
    return 0;
  }
  // deform
  int deform() {
    shared_ptr<deform_energy> e = std::make_shared<deform_energy>(tris_, nods_, cross_);

    Map<Matrix<double, 6, 1>> X(nods_.begin());
    QUERY_ENERGY_VALUE(e, X);

    VectorXd rhs = VectorXd::Zero(6);
    e->gra(&X[0], rhs.data());
    rhs = -rhs;
    Matrix<double, 6, 6> H;
    H.setZero();
    e->hes(nullptr, H.data());
    SparseMatrix<double> Hs = H.sparseView();

    vector<size_t> g2l{static_cast<size_t>(-1), static_cast<size_t>(-1), 0, 1, 2, 3};
    surfparam::rm_spmat_col_row(Hs, g2l);
    surfparam::rm_vector_row(rhs, g2l);
    SimplicialCholesky<SparseMatrix<double>> sol;
    sol.compute(Hs);
    ASSERT(sol.info() == Success);
    VectorXd dx = sol.solve(rhs);
    ASSERT(sol.info() == Success);
    VectorXd Dx = VectorXd::Zero(6);
    surfparam::rc_vector_row(dx, g2l, Dx);
    X += Dx;
    QUERY_ENERGY_VALUE(e, X);

    return 0;
  }
  // debug
  int see_triangle_height(const char *filename) const {
    matd_t ht(2, 3);
    calc_barycentric_grad(nods_.begin(), ht.begin());
    mati_t line(2, 3);
    line(0, colon()) = colon(0, 2);
    line(1, colon()) = colon(3, 5);
    matd_t vert(2, 6);
    vert(colon(), colon(0, 2)) = nods_;
    vert(colon(), colon(3, 5)) = vert(colon(), colon(0, 2))-ht;
    matd_t v = zeros<double>(3, 6);
    v(colon(0, 1), colon()) = vert;
    ofstream os(filename);
    line2vtk(os, v.begin(), v.size(2), line.begin(), line.size(2));
    os.close();
    return 0;
  }
  int run_test() const {
    matd_t ht(2, 3);
    calc_barycentric_grad(nods_.begin(), ht.begin());
    cout << ht*ones<double>(3, 1) << endl;
    return 0;
  }
  int see_scalar_filed_grad(const char *filename) {
    return 0;
  }
private:
  mati_t tris_;    // 3x1
  matd_t nods_;    // 2x3
  Matrix2d cross_; // 2x2
};

int main(int argc, char *argv[])
{
  boost::filesystem::create_directory("./cross_param");

  cross_field_param cfp;
  cfp.set_triangle();
  cfp.set_cross();
  cfp.see_triangle("./cross_param/origin.vtk");
  cfp.see_cross("./cross_param/cross.vtk");
//  cfp.run_test();

  cfp.deform();
  cfp.see_triangle("./cross_param/deform.vtk");

  cout << "[info] done\n";
  return 0;
}
