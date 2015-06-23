#include "deform_transfer.h"

#include <jtflib/mesh/io.h>
#include <jtflib/mesh/mesh.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <zjucad/matrix/io.h>
#include <Eigen/Geometry>

#include "util.h"
#include "vtk.h"
#include "nanoflann.hpp"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;
using namespace surfparam;
using namespace jtf::mesh;
using namespace nanoflann;

namespace geom_deform {

extern "C" {

void unit_deform_energy_(double *val, const double *x, const double *Tinv, const double *S);
void unit_deform_energy_jac_(double *jac, const double *x, const double *Tinv, const double *S);
void unit_deform_energy_hes_(double *hes, const double *x, const double *Tinv, const double *S);

void unit_smooth_energy_(double *val, const double *x, const double *di, const double *dj);
void unit_smooth_energy_jac_(double *jac, const double *x, const double *di, const double *dj);
void unit_smooth_energy_hes_(double *hes, const double *x, const double *di, const double *dj);

void unit_identity_energy_(double *val, const double *x, const double *d);
void unit_identity_energy_jac_(double *jac, const double *x, const double *d);
void unit_identity_energy_hes_(double *hes, const double *x, const double *d);

}

int deform_transfer::debug_unit_energy() const {
  srand(time(NULL));
  matd_t nods = rand<double>(3, 4);
  matd_t d = nods(colon(), colon(1, 3))-nods(colon(), 0)*ones<double>(1, 3);
  matd_t temp = d;
  if ( inv(d) ) {
    cerr << "[info] not invertible\n";
    return __LINE__;
  }
  cout << temp*d << endl;

  double value = 0;
  unit_identity_energy_(&value, &nods[0], &d[0]);
  cout << "value: " << value << endl;

  matd_t grad = zeros<double>(12, 1);
  unit_identity_energy_jac_(&grad[0], &nods[0], &d[0]);
  cout << "grad: " << grad << endl;

  matd_t H(12, 12);
  unit_identity_energy_hes_(&H[0], &nods[0], &d[0]);
  cout << "hes: " << H << endl;

  return 0;
}

class dt_deform_energy : public Functional<double>
{
public:
  typedef zjucad::matrix::matrix<size_t> mati_t;
  typedef zjucad::matrix::matrix<double> matd_t;
  dt_deform_energy(const mati_t &tar_cell, const matd_t &tar_nods,
                   const MatrixXd &Sinv, const MatrixXd &Tinv,
                   const set<tuple<size_t, size_t>> mapping, double w)
    : tris_(tar_cell),
      nods_(tar_nods),
      Sinv_(Sinv),
      Tinv_(Tinv),
      mapping_(mapping),
      w_(w) {}
  size_t Nx() const {
    return nods_.size();
  }
  int Val(const double *x, double *val) const {
    itr_matrix<const double *> X(3, Nx()/3, x);
    for (auto &e : mapping_) {
      const size_t src_fa = std::get<0>(e);
      const size_t tar_fa = std::get<1>(e);
      matd_t vert = X(colon(), tris_(colon(), tar_fa));
      double value = 0;
      unit_deform_energy_(&value, &vert[0], &Tinv_(0, 3*tar_fa), &src_grad_(0, 3*src_fa));
      *val += w_*value;
    }
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    itr_matrix<const double *> X(3, Nx()/3, x);
    itr_matrix<double *> G(3, Nx()/3, gra);
    for (auto &e : mapping_) {
      const size_t src_fa = std::get<0>(e);
      const size_t tar_fa = std::get<1>(e);
      matd_t vert = X(colon(), tris_(colon(), tar_fa));
      matd_t g = zeros<double>(3, 4);
      unit_deform_energy_jac_(&g[0], &vert[0], &Tinv_(0, 3*tar_fa), &src_grad_(0, 3*src_fa));
      for (size_t j = 0; j < 4; ++j)
        G(colon(), tris_(j, tar_fa)) += w_*g(colon(), j);
    }
    return 0;
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    for (auto &e : mapping_) {
      const size_t src_fa = std::get<0>(e);
      const size_t tar_fa = std::get<1>(e);
      matd_t H = zeros<double>(12, 12);
      unit_deform_energy_hes_(&H[0], NULL, &Tinv_(0, 3*tar_fa), &src_grad_(0, 3*src_fa));
      for (size_t p = 0; p < 12; ++p) {
        for (size_t q = 0; q < 12; ++q) {
          const size_t I = 3*tris_(p/3, tar_fa)+p%3;
          const size_t J = 3*tris_(q/3, tar_fa)+q%3;
          if ( H(p, q) != 0.0 )
            hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
        }
      }
    }
    return 0;
  }
  void UpdateSourceGrad(const mati_t &src_tris, const matd_t &src_def) {
    src_grad_.resize(3, 3*src_tris.size(2));
#pragma omp parallel for
    for (size_t i = 0; i < src_tris.size(2); ++i) {
      matd_t base = src_def(colon(), src_tris(colon(1, 3), i))
          - src_def(colon(), src_tris(0, i))*ones<double>(1, 3);
      src_grad_.block<3, 3>(0, 3*i) = Map<const Matrix3d>(&base[0])*Sinv_.block<3, 3>(0, 3*i);
    }
  }
private:
  const mati_t &tris_;
  const matd_t &nods_;
  const MatrixXd &Sinv_;
  const MatrixXd &Tinv_;
  const set<tuple<size_t, size_t>> &mapping_;
  double w_;
  MatrixXd src_grad_;
};

class dt_smooth_energy : public Functional<double>
{
public:
  typedef zjucad::matrix::matrix<size_t> mati_t;
  typedef zjucad::matrix::matrix<double> matd_t;
  dt_smooth_energy(const mati_t &src_cell, const matd_t &src_nods,
                   const MatrixXd &Sinv, const double w)
    : tris_(src_cell), nods_(src_nods), Sinv_(Sinv), w_(w) {
    mati_t tris = tris_(colon(0, 2), colon());
    e2c_.reset(edge2cell_adjacent::create(tris, false));
  }
  size_t Nx() const {
    return nods_.size();
  }
  int Val(const double *x, double *val) const {
    itr_matrix<const double *> X(3, Nx()/3, x);
    for (auto &e : e2c_->edges_) {
      pair<size_t, size_t> face = e2c_->query(e.first, e.second);
      if ( e2c_->is_boundary_edge(face) )
        continue;
      const size_t fa[2] = {face.first, face.second};
      matd_t vert(3, 8);
      vert(colon(), colon(0, 3)) = X(colon(), tris_(colon(), fa[0]));
      vert(colon(), colon(4, 7)) = X(colon(), tris_(colon(), fa[1]));
      double value = 0;
      unit_smooth_energy_(&value, &vert[0], &Sinv_(0, 3*fa[0]), &Sinv_(0, 3*fa[1]));
      *val += w_* value;
    }
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    itr_matrix<const double *> X(3, Nx()/3, x);
    itr_matrix<double *> grad(3, Nx()/3, gra);
    for (auto &e : e2c_->edges_) {
      pair<size_t, size_t> face = e2c_->query(e.first, e.second);
      if ( e2c_->is_boundary_edge(face) )
        continue;
      const size_t fa[2] = {face.first, face.second};
      matd_t vert(3, 8);
      vert(colon(), colon(0, 3)) = X(colon(), tris_(colon(), fa[0]));
      vert(colon(), colon(4, 7)) = X(colon(), tris_(colon(), fa[1]));
      matd_t g = zeros<double>(3, 8);
      unit_smooth_energy_jac_(&g[0], &vert[0], &Sinv_(0, 3*fa[0]), &Sinv_(0, 3*fa[1]));
      for (size_t j = 0; j < 8; ++j) {
        grad(colon(), tris_(j%4, fa[j/4])) += w_*g(colon(), j);
      }
    }
    return 0;
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    for (auto &e : e2c_->edges_) {
      pair<size_t, size_t> face = e2c_->query(e.first, e.second);
      if ( e2c_->is_boundary_edge(face) )
        continue;
      const size_t fa[2] = {face.first, face.second};
      matd_t H = zeros<double>(24, 24);
      unit_smooth_energy_hes_(&H[0], NULL, &Sinv_(0, 3*fa[0]), &Sinv_(0, 3*fa[1]));
      for (size_t p = 0; p < 24; ++p) {
        for (size_t q = 0; q < 24; ++q) {
          const size_t I = 3*tris_((p/3)%4, fa[p/12])+p%3;
          const size_t J = 3*tris_((q/3)%4, fa[q/12])+q%3;
          if ( H(p, q) != 0.0 )
            hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
        }
      }
    }
    return 0;
  }
  void ResetWeight(const double w) {
    w_ = w;
  }
private:
  const mati_t &tris_;
  const matd_t &nods_;
  double w_;
  shared_ptr<edge2cell_adjacent> e2c_;
  const MatrixXd &Sinv_;
};

class dt_identity_energy : public Functional<double>
{
public:
  typedef zjucad::matrix::matrix<size_t> mati_t;
  typedef zjucad::matrix::matrix<double> matd_t;
  dt_identity_energy(const mati_t &src_cell, const matd_t &src_nods,
                     const MatrixXd &Sinv, const double w)
    : tris_(src_cell), nods_(src_nods), Sinv_(Sinv), w_(w) {}
  size_t Nx() const {
    return nods_.size();
  }
  int Val(const double *x, double *val) const {
    itr_matrix<const double *> X(3, Nx()/3, x);
    for (size_t i = 0; i < tris_.size(2); ++i) {
      matd_t vert = X(colon(), tris_(colon(), i));
      double value = 0;
      unit_identity_energy_(&value, &vert[0], &Sinv_(0, 3*i));
      *val += w_*value;
    }
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    itr_matrix<const double *> X(3, Nx()/3, x);
    itr_matrix<double *> grad(3, Nx()/3, gra);
    for (size_t i = 0; i < tris_.size(2); ++i) {
      matd_t vert = X(colon(), tris_(colon(), i));
      matd_t g = zeros<double>(3, 4);
      unit_identity_energy_jac_(&g[0], &vert[0], &Sinv_(0, 3*i));
      for (size_t j = 0; j < 4; ++j) {
        grad(colon(), tris_(j, i)) += w_*g(colon(), j);
      }
    }
    return 0;
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    for (size_t i = 0; i < tris_.size(2); ++i) {
      matd_t H = zeros<double>(12, 12);
      unit_identity_energy_hes_(&H[0], NULL, &Sinv_(0, 3*i));
      for (size_t p = 0; p < 12; ++p) {
        for (size_t q = 0; q < 12; ++q) {
          const size_t I = 3*tris_(p/3, i)+p%3;
          const size_t J = 3*tris_(q/3, i)+q%3;
          if ( H(p, q) != 0.0 )
            hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
        }
      }
    }
    return 0;
  }
  void ResetWeight(const double w) {
    w_ = w;
  }
private:
  const mati_t &tris_;
  const matd_t &nods_;
  double w_;
  const MatrixXd &Sinv_;
};

class dt_distance_energy : public Functional<double>
{
public:
  typedef zjucad::matrix::matrix<size_t> mati_t;
  typedef zjucad::matrix::matrix<double> matd_t;
  typedef KDTreeEigenMatrixAdaptor<Matrix<double, -1, -1>> kd_tree_t;
  dt_distance_energy(const mati_t &src_cell, const matd_t &src_nods,
                     const matd_t &tar_nods, const double w)
    : tris_(src_cell), nods_(src_nods), w_(w) {
    Map<const MatrixXd> temp(&tar_nods[0], tar_nods.size(1), tar_nods.size(2));
    pts_ = temp.transpose();
    kdt_ = std::make_shared<kd_tree_t>(3, pts_, 10);
    kdt_->index->buildIndex();
    c_ = zeros<double>(src_nods.size(1), src_nods.size(2));
  }
  size_t Nx() const {
    return nods_.size();
  }
  int Val(const double *x, double *val) const {
    itr_matrix<const double *> X(3, Nx()/3, x);
    for (size_t i = 0; i < X.size(2); ++i) {
      matd_t diff = X(colon(), i) - c_(colon(), i);
      *val += w_*dot(diff, diff);
    }
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    itr_matrix<const double *> X(3, Nx()/3, x);
    itr_matrix<double *> G(3, Nx()/3, gra);
#pragma omp parallel for
    for (size_t i = 0; i < X.size(2); ++i) {
      G(colon(), i) += 2*w_*(X(colon(), i)-c_(colon(), i));
    }
    return 0;
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    for (size_t i = 0; i < Nx(); ++i)
      hes->push_back(Triplet<double>(i, i, 2*w_));
    return 0;
  }
  void ResetWeight(const double w) {
    w_ = w;
  }
  int UpdateClosetPoints(const double *x) {
    itr_matrix<const double *> X(3, Nx()/3, x);
    const size_t num_results = 3;
#pragma omp parallel for
    for (size_t i = 0; i < X.size(2); ++i) {
      vector<size_t> ret_idx(num_results);
      vector<double> sqr_dist(num_results);
      KNNResultSet<double> result_set(num_results);
      result_set.init(&ret_idx[0], &sqr_dist[0]);
      kdt_->index->findNeighbors(result_set, &X(0, i), SearchParams(10));
      ASSERT(ret_idx[0] >= 0 && ret_idx[0] < pts_.rows());
      c_(0, i) = pts_(ret_idx[0], 0);
      c_(1, i) = pts_(ret_idx[0], 1);
      c_(2, i) = pts_(ret_idx[0], 2);
    }
    return 0;
  }
private:
  bool is_valid() const;

  const mati_t &tris_;
  const matd_t &nods_;
  MatrixXd pts_;
  shared_ptr<kd_tree_t> kdt_;
  matd_t c_;
  double w_;
};

int deform_transfer::debug_energies() const {
  cout << "[info] debug energy functional\n";
  shared_ptr<Functional<double>> e0
      = std::make_shared<dt_smooth_energy>(src_tris_, src_ref_nods_, Sinv_, 1.0);
  shared_ptr<Functional<double>> e1
      = std::make_shared<dt_identity_energy>(src_tris_, src_ref_nods_, Sinv_, 1.0);
  {
    matd_t x = src_ref_nods_;
    double val = 0;
    e0->Val(&x[0], &val);
    cout << "energy value: " << val << endl;
    matd_t grad = zeros<double>(e0->Nx(), 1);
    e0->Gra(&x[0], &grad[0]);
    cout << "grad max: " << max(grad) << endl << endl;
  }
  {
    matd_t x = src_ref_nods_;
    Matrix3d rot;
    rot = AngleAxisd(M_PI/2, Vector3d::UnitX());
    matd_t R(3, 3);
    std::copy(rot.data(), rot.data()+9, R.begin());
    x = temp(R*x);
    double val = 0;
    e0->Val(&x[0], &val);
    cout << "energy value: " << val << endl;
    matd_t grad = zeros<double>(e0->Nx(), 1);
    e0->Gra(&x[0], &grad[0]);
    cout << "grad max: " << max(grad) << endl << endl;
  }
  {
    matd_t x = src_ref_nods_;
    double val = 0;
    e1->Val(&x[0], &val);
    cout << "energy value: " << val << endl;
    matd_t grad = zeros<double>(e1->Nx(), 1);
    e1->Gra(&x[0], &grad[0]);
    cout << "grad max: " << max(grad) << endl << endl;
  }
  {
    matd_t x = 2*src_ref_nods_;
    double val = 0;
    e1->Val(&x[0], &val);
    cout << "energy value: " << val << endl;
    matd_t grad = zeros<double>(e1->Nx(), 1);
    e1->Gra(&x[0], &grad[0]);
    cout << "grad max: " << norm(grad) << endl << endl;
  }

  return 0;
}

void deform_transfer::append_fourth_vert(const mati_t &tri_cell, const matd_t &tri_nods,
                                         mati_t &tet_cell, matd_t &tet_nods) const {
  const size_t nbr_tri_vert = tri_nods.size(2);
  const size_t nbr_tet_vert = nbr_tri_vert + tri_cell.size(2);
  tet_cell.resize(4, tri_cell.size(2));
  tet_cell(colon(0, 2), colon()) = tri_cell;
  tet_cell(3, colon()) = colon(nbr_tri_vert, nbr_tet_vert-1);
  tet_nods.resize(3, nbr_tet_vert);
  tet_nods(colon(), colon(0, nbr_tri_vert-1)) = tri_nods;
#pragma omp parallel for
  for (size_t i = 0; i < tri_cell.size(2); ++i) {
    matd_t vert = tri_nods(colon(), tri_cell(colon(), i));
    matd_t temp = cross(vert(colon(), 1)-vert(colon(), 0), vert(colon(), 2)-vert(colon(), 0));
    tet_nods(colon(), nbr_tri_vert+i) = vert(colon(), 0) + temp/std::sqrt(norm(temp));
  }
}

void deform_transfer::remove_fourth_vert(const mati_t &tet_cell, const matd_t &tet_nods,
                                         mati_t &tri_cell, matd_t &tri_nods) const {
  tri_cell.resize(3, tet_cell.size(2));
  tri_cell = tet_cell(colon(0, 2), colon());
  const size_t nbr_tri_vert = max(tri_cell)+1;
  tri_nods.resize(3, nbr_tri_vert);
  tri_nods = tet_nods(colon(), colon(0, nbr_tri_vert-1));
}

deform_transfer::deform_transfer() {}

int deform_transfer::load_reference_source_mesh(const char *filename) {
  mati_t tris;
  matd_t nods;
  int rtn = jtf::mesh::load_obj(filename, tris, nods);
  append_fourth_vert(tris, nods, src_tris_, src_ref_nods_);
  return rtn;
}

int deform_transfer::load_reference_target_mesh(const char *filename) {
  mati_t tris;
  matd_t nods;
  int rtn = jtf::mesh::load_obj(filename, tris, nods);
  append_fourth_vert(tris, nods, tar_tris_, tar_ref_nods_);
  return rtn;
}

/// can be invoked multiple times
int deform_transfer::load_deformed_source_mesh(const char *filename) {
  mati_t tris;
  matd_t nods;
  int rtn = jtf::mesh::load_obj(filename, tris, nods);
  append_fourth_vert(tris, nods, src_tris_, src_def_nods_);
  return rtn;
}

int deform_transfer::load_vertex_markers(const char *filename) {
  FILE *fp = fopen(filename, "r");
  if ( fp == NULL ) {
    cerr << "[error] can not open markers\n";
    return __LINE__;
  }
  size_t nbr_markers;
  int state = fscanf(fp, "%zu", &nbr_markers);
  cout << "[info] number of markers: " << nbr_markers << endl;
  vert_map_.resize(nbr_markers);
  for (size_t i = 0; i < nbr_markers; ++i) {
    size_t src_mark, tar_mark;
    state = fscanf(fp, "%zu, %zu", &src_mark, &tar_mark);
    vert_map_[i] = std::make_tuple(src_mark, tar_mark);
  }
  fclose(fp);
  return 0;
}

int deform_transfer::see_target_markers(const char *filename) const {
  const size_t nbr_marker = vert_map_.size();
  mati_t point = colon(0, nbr_marker-1);
  matd_t nods(3, nbr_marker);
#pragma omp parallel for
  for (size_t i = 0; i < vert_map_.size(); ++i)
    nods(colon(), i) = tar_ref_nods_(colon(), std::get<1>(vert_map_[i]));
  ofstream os(filename);
  point2vtk(os, nods.begin(), nods.size(2), point.begin(), point.size());
  return 0;
}

int deform_transfer::save_reference_source_mesh(const char *filename) const {
  mati_t tris;
  matd_t nods;
  remove_fourth_vert(src_tris_, src_ref_nods_, tris, nods);
  return jtf::mesh::save_obj(filename, tris, nods);
}

int deform_transfer::save_reference_target_mesh(const char *filename) const {
  mati_t tris;
  matd_t nods;
  remove_fourth_vert(tar_tris_, tar_ref_nods_, tris, nods);
  return jtf::mesh::save_obj(filename, tris, nods);
}

int deform_transfer::save_deformed_source_mesh(const char *filename) const {
  mati_t tris;
  matd_t nods;
  remove_fourth_vert(src_tris_, src_def_nods_, tris, nods);
  return jtf::mesh::save_obj(filename, tris, nods);
}

int deform_transfer::save_deformed_target_mesh(const char *filename) const {
  mati_t tris;
  matd_t nods;
  remove_fourth_vert(tar_tris_, tar_def_nods_, tris, nods);
  return jtf::mesh::save_obj(filename, tris, nods);
}

int deform_transfer::see_ghost_tet_mesh(const char *filename, const string &which) const {
  ofstream os(filename);
  if ( which == "source_ref" )
    tet2vtk(os, &src_ref_nods_[0], src_ref_nods_.size(2), &src_tris_[0], src_tris_.size(2));
  else if ( which == "target_ref" )
    tet2vtk(os, &tar_ref_nods_[0], tar_ref_nods_.size(2), &tar_tris_[0], tar_tris_.size(2));
  else
    return __LINE__;
  return 0;
}

int deform_transfer::see_corres_mesh(const char *filename) const {
  mati_t tris;
  matd_t nods;
  remove_fourth_vert(src_tris_, src_cor_nods_, tris, nods);
  return jtf::mesh::save_obj(filename, tris, nods);
}

int deform_transfer::init() {
  cout << "[info] initialization\n";
  // calculate Sinv
  Sinv_.resize(3, 3*src_tris_.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < src_tris_.size(2); ++i) {
    matd_t base = src_ref_nods_(colon(), src_tris_(colon(1, 3), i))
        - src_ref_nods_(colon(), src_tris_(0, i))*ones<double>(1, 3);
    Matrix3d LB(&base[0]);
    FullPivLU<Matrix3d> lu(LB);
    if ( lu.isInvertible() ) {
      Sinv_.block<3, 3>(0, 3*i) = lu.inverse();
    } else {
      cerr << "[error] degenerated triagnle in source mesh\n";
      exit(EXIT_FAILURE);
    }
  }
  // calculate Tinv
  Tinv_.resize(3, 3*tar_tris_.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < tar_tris_.size(2); ++i) {
    matd_t base = tar_ref_nods_(colon(), tar_tris_(colon(1, 3), i))
        - tar_ref_nods_(colon(), tar_tris_(0, i))*ones<double>(1, 3);
    Matrix3d LB(&base[0]);
    FullPivLU<Matrix3d> lu(LB);
    if ( lu.isInvertible() ) {
      Tinv_.block<3, 3>(0, 3*i) = lu.inverse();
    } else {
      cerr << "[error] degenerated triangle in target mesh\n";
      exit(EXIT_FAILURE);
    }
  }
  return 0;
}

int deform_transfer::solve_corres_precompute() {
  cout << "[info] precomputation for correspondence\n";
  // handle constraints
  for (auto &e : vert_map_) {
    const size_t id = std::get<0>(e);
    fix_dof_.insert(3*id+0);
    fix_dof_.insert(3*id+1);
    fix_dof_.insert(3*id+2);
  }
  g2l_.resize(src_ref_nods_.size());
  size_t ptr = 0;
  for (size_t i = 0; i < g2l_.size(); ++i) {
    if ( fix_dof_.find(i) != fix_dof_.end() )
      g2l_[i] = -1;
    else
      g2l_[i] = ptr++;
  }
  return 0;
}

int deform_transfer::solve_corres_first_phase() {
  cout << "[info] first phase for resolving correspondence\n";
  // assemble energy
  vector<double> w{2.0, 0.001, 0.0};
  buff_.resize(3);
  buff_[SMOOTH] = std::make_shared<dt_smooth_energy>(src_tris_, src_ref_nods_, Sinv_, w[SMOOTH]);
  buff_[IDENTITY] = std::make_shared<dt_identity_energy>(src_tris_, src_ref_nods_, Sinv_, w[IDENTITY]);
  buff_[DISTANCE] = std::make_shared<dt_distance_energy>(src_tris_, src_ref_nods_, tar_ref_nods_, w[DISTANCE]);
  try {
    corre_e_ = std::make_shared<energy_t<double>>(buff_);
  } catch ( exception &e ) {
    cerr << "[exception] " << e.what() << endl;
    exit(EXIT_FAILURE);
  }

  /// set initial value
  cout << "\t@set initial values\n";
  src_cor_nods_ = src_ref_nods_;
#pragma omp parallel for
  for (size_t i = 0; i < vert_map_.size(); ++i) {
    src_cor_nods_(colon(), std::get<0>(vert_map_[i]))
        = tar_ref_nods_(colon(), std::get<1>(vert_map_[i]));
  }

  /// solve
  const size_t dim = corre_e_->Nx();
  Map<VectorXd> x(&src_cor_nods_[0], dim);

  double energy_prev = 0;
  corre_e_->Val(&x[0], &energy_prev);
  cout << "[info] prev energy value: " << energy_prev << endl;

  vector<Triplet<double>> trips;
  corre_e_->Hes(&x[0], &trips);
  SparseMatrix<double> H(dim, dim);
  H.reserve(trips.size());
  H.setFromTriplets(trips.begin(), trips.end());

  VectorXd rhs = VectorXd::Zero(dim);
  corre_e_->Gra(&x[0], rhs.data());
  rhs = -rhs;

  if ( !fix_dof_.empty() ) {
    rm_spmat_col_row(H, g2l_);
    rm_vector_row(rhs, g2l_);
  }

  SimplicialCholesky<SparseMatrix<double>> sol;
  sol.compute(H);
  ASSERT(sol.info() == Success);
  VectorXd dx = sol.solve(rhs);
  ASSERT(sol.info() == Success);

  VectorXd Dx = VectorXd::Zero(dim);
  if ( !fix_dof_.empty() ) {
    rc_vector_row(dx, g2l_, Dx);
  } else {
    Dx = dx;
  }
  x += Dx;

  double energy_post = 0;
  corre_e_->Val(&x[0], &energy_post);
  cout << "[info] post energy value: " << energy_post << endl;

  cout << "[info] first phase completed\n";
  return 0;
}

int deform_transfer::solve_corres_second_phase() {
  cout << "[info] second phase for resolving correspondence\n";

  const size_t dim = corre_e_->Nx();
  Map<VectorXd> x(&src_cor_nods_[0], dim);
  SimplicialCholesky<SparseMatrix<double>> sol;

  for (size_t iter = 0; iter < 4; ++iter) {
    cout << "[info] iter " << iter << endl;
    const double wc = 1.0 + iter*(5000.0-1.0)/3;
    cout << "[info] current weight: " << wc << endl;
    std::dynamic_pointer_cast<dt_distance_energy>(buff_[DISTANCE])
        ->ResetWeight(wc);
    std::dynamic_pointer_cast<dt_distance_energy>(buff_[DISTANCE])
        ->UpdateClosetPoints(&x[0]);

    vector<Triplet<double>> trips;
    corre_e_->Hes(&x[0], &trips);
    SparseMatrix<double> H(dim, dim);
    H.reserve(trips.size());
    H.setFromTriplets(trips.begin(), trips.end());

    VectorXd rhs = VectorXd::Zero(dim);
    corre_e_->Gra(&x[0], rhs.data());
    rhs = -rhs;

    if ( !fix_dof_.empty() ) {
      rm_spmat_col_row(H, g2l_);
      rm_vector_row(rhs, g2l_);
    }

    sol.compute(H);
    ASSERT(sol.info() == Success);
    VectorXd dx = sol.solve(rhs);
    ASSERT(sol.info() == Success);

    VectorXd Dx = VectorXd::Zero(dim);
    if ( !fix_dof_.empty() ) {
      rc_vector_row(dx, g2l_, Dx);
    } else {
      Dx = dx;
    }
    x += Dx;
  }
  cout << "[info] second phase completed\n";
  return 0;
}

void deform_transfer::build_corre_face(const mati_t &src_tris, const matd_t &src_nods,
                                       const mati_t &tar_tris, const matd_t &tar_nods,
                                       set<tuple<size_t, size_t>> &mappings, bool need_swap) {}

int deform_transfer::compute_triangle_corres() {
  cout << "[info] computing triangle correspondence...";
  typedef KDTreeEigenMatrixAdaptor<MatrixXd> kd_tree_t;
  matd_t src_cent(3, src_tris_.size(2)), tar_cent(3, tar_tris_.size(2));
#pragma omp parallel for
  for (size_t i = 0; i < src_cent.size(2); ++i) {
    src_cent(colon(), i) = src_cor_nods_(colon(), src_tris_(colon(0, 2), i))*ones<double>(3, 1)/3.0;
  }
#pragma omp parallel for
  for (size_t i = 0; i < tar_cent.size(2); ++i) {
    tar_cent(colon(), i) = tar_ref_nods_(colon(), tar_tris_(colon(0, 2), i))*ones<double>(3, 1)/3.0;
  }
  {
    MatrixXd pts = Map<const MatrixXd>(&tar_cent[0], tar_cent.size(1), tar_cent.size(2)).transpose();
    kd_tree_t kdt(3, pts, 10);
    kdt.index->buildIndex();
#pragma omp parallel for
    for (size_t i = 0; i < src_cent.size(2); ++i) {
      // do a radius search
    }
  }
  {

  }
  cout << "...complete\n";
  return 0;
}

int deform_transfer::deformation_transfer_precompute() {
  cout << "[info] precomputation for deformation transfer...";
  tar_def_nods_ = tar_ref_nods_;
  deform_e_ = std::make_shared<dt_deform_energy>(tar_tris_, tar_ref_nods_, Sinv_, Tinv_, tri_map_, 1.0);
  cout << "complete\n";
  return 0;
}

int deform_transfer::deformation_transfer() {
  cout << "[info] deformation transfer begin...";
  std::dynamic_pointer_cast<dt_deform_energy>(deform_e_)
      ->UpdateSourceGrad(src_tris_, src_def_nods_);

  const size_t dim = deform_e_->Nx();
  Map<VectorXd> x(&tar_def_nods_[0], dim);

  vector<Triplet<double>> trips;
  deform_e_->Hes(&x[0], &trips);
  SparseMatrix<double> H(dim, dim);
  H.reserve(trips.size());
  H.setFromTriplets(trips.begin(), trips.end());

  VectorXd rhs = VectorXd::Zero(dim);
  deform_e_->Gra(&x[0], rhs.data());
  rhs = -rhs;

  SimplicialCholesky<SparseMatrix<double>> sol;
  sol.compute(H);
  ASSERT(sol.info() == Success);
  VectorXd dx = sol.solve(rhs);
  ASSERT(sol.info() == Success);

  x += dx;

  cout << "complete\n";
  return 0;
}

}
