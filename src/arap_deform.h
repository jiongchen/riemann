#ifndef ARAP_DEFORM_H
#define ARAP_DEFORM_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <zjucad/matrix/matrix.h>
#include <unordered_set>

namespace core {

class arap_energy;

class arap_deform
{
public:
  typedef zjucad::matrix::matrix<size_t> mati_t;
  typedef zjucad::matrix::matrix<double> matd_t;
  arap_deform(const mati_t &tris, const matd_t &nods);
  int pre_compute(const std::vector<size_t> &idx);
  int deformation(double *x);
private:
  const mati_t tris_;
  const matd_t nods_;

  std::unordered_set<size_t> fixed_dofs_;
  std::vector<size_t> g2l_;

  std::shared_ptr<arap_energy> e_;
  Eigen::SparseMatrix<double> L_;
  Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> sol_;
};

}

#endif
