#include "deform_transfer.h"

#include <jtflib/mesh/io.h>
#include <jtflib/mesh/mesh.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <zjucad/matrix/io.h>
#include <Eigen/Geometry>

#include "util.h"
#include "vtk.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;
using namespace surfparam;
using namespace jtf::mesh;

namespace geom_deform {

extern "C" {

void unit_deform_energy_();
void unit_deform_energy_jac_();
void unit_deform_energy_hes_();

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
  dt_deform_energy(const mati_t &src_cell, const matd_t &src_nods, const double w);
  size_t Nx() const {

  }
  int Val(const double *x, double *val) const {

  }
  int Gra(const double *x, double *gra) const {

  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {

  }
  void ResetWeight(const double w) {
    w_ = w;
  }
private:
  const mati_t &tris_;
  const matd_t &nods_;
  double w_;
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
      const size_t fa[2] = {e.first, e.second};
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
      const size_t fa[2] = {e.first, e.second};
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
      const size_t fa[2] = {e.first, e.second};
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
  dt_distance_energy(const mati_t &src_cell, const matd_t &src_nods, const double w)
    : tris_(src_cell), nods_(src_nods), w_(w) {}
  size_t Nx() const {
    return nods_.size();
  }
  int Val(const double *x, double *val) const {
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    return 0;
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    return 0;
  }
  void ResetWeight(const double w) {
    w_ = w;
  }
private:
  const mati_t &tris_;
  const matd_t &nods_;
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

int deform_transfer::solve_corres_precompute() {
  cout << "[info] precomputation for correspondence\n";
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
      cerr << "[error] degenerated triagnle\n";
      exit(EXIT_FAILURE);
    }
  }
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
  buff_[DISTANCE] = std::make_shared<dt_distance_energy>(src_tris_, src_ref_nods_, w[DISTANCE]);
  try {
    corre_e_ = std::make_shared<energy_t<double>>(buff_);
  } catch ( exception &e ) {
    cerr << "[exception] " << e.what() << endl;
    exit(EXIT_FAILURE);
  }

  // set initial value
  cout << "\t@set initial values\n";
  src_cor_nods_ = src_ref_nods_;
#pragma omp parallel for
  for (size_t i = 0; i < vert_map_.size(); ++i) {
    src_cor_nods_(colon(), std::get<0>(vert_map_[i]))
        = tar_ref_nods_(colon(), std::get<1>(vert_map_[i]));
  }

  // solve
  cout << "\t@solve\n";
  const size_t dim = corre_e_->Nx();
  ASSERT(src_cor_nods_.size() == dim);
  Map<VectorXd> x(&src_cor_nods_[0], dim);
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

  cout << "[info] first phase completed\n";
  return 0;
}

}
