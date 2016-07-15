#include "volume_frame.h"

#include <zjucad/matrix/itr_matrix.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <zjucad/linear_solver/linear_solver.h>
#include <Eigen/Eigen>

#ifdef USE_CHOLMOD
  #include <Eigen/CholmodSupport>
#endif

#include "config.h"
#include "grad_operator.h"
#include "lbfgs_solve.h"
#include "sh_zyz_convert.h"
#include "util.h"

using namespace std;
using namespace zjucad::matrix;
using namespace jtf::mesh;
using namespace Eigen;

namespace riemann {

extern "C" {

  void cubic_sym_smooth_(double *val, const double *abc, const double *CR, const double *vol);
  void cubic_sym_smooth_jac_(double *jac, const double *abc, const double *CR, const double *vol);
    
  void cubic_sym_align_(double *val, const double *abc, const double *rnz, const double *area);
  void cubic_sym_align_jac_(double *jac, const double *abc, const double *rnz, const double *area);

  void cubic_smooth_sh_coef_(double *val, const double *F, const double *CR, const double *vol);
  void cubic_smooth_sh_coef_jac_(double *jac, const double *F, const double *CR, const double *vol);
  void cubic_smooth_sh_coef_hes_(double *hes, const double *F, const double *CR, const double *vol);

  void cubic_align_sh_coef_(double *val, const double *F, const double *Rnz, const double *area);
  void cubic_align_sh_coef_jac_(double *jac, const double *F, const double *Rnz, const double *area);
  void cubic_align_sh_coef_hes_(double *hes, const double *F, const double *Rnz, const double *area);

  void cubic_sym_smooth_tet_(double *val, const double *abc, const double *stiff);
  void cubic_sym_smooth_tet_jac_(double *jac, const double *abc, const double *stiff);
  
}

static inline void normal2zyz(const double *n, double *zyz);

//===============================================================================

class SH_smooth_energy : public Functional<double>
{
public:
  SH_smooth_energy(const mati_t &tets, const matd_t &nods, const double w);
  size_t Nx() const;
  int Val(const double *abc, double *val) const;
  int Gra(const double *abc, double *gra) const;
  int Hes(const double *abc, std::vector<Eigen::Triplet<double>> *hes) const { return __LINE__; }
  int ValSH(const double *f, double *val) const;
  int GraSH(const double *f, double *gra) const;
  int HesSH(const double *f, std::vector<Eigen::Triplet<double>> *hes) const;
private:
  double w_;
  const size_t dim_;
  const mati_t &tets_;
  matd_t basis_;
  matd_t vol_;
};

class SH_align_energy : public Functional<double>
{
public:
  SH_align_energy(const mati_t &tets, const matd_t &nods, const double w);
  size_t Nx() const;
  int Val(const double *abc, double *val) const;
  int Gra(const double *abc, double *gra) const;
  int Hes(const double *abc, std::vector<Eigen::Triplet<double>> *hes) const { return __LINE__; }
  int ValSH(const double *f, double *val) const;
  int GraSH(const double *f, double *gra) const;
  int HesSH(const double *f, std::vector<Eigen::Triplet<double>> *hes) const;
private:
  double w_;
  const size_t dim_;
  const mati_t &tets_;
  mati_t surf_;
  matd_t zyz_;
  matd_t area_;
};

class SH_smooth_energy_tet : public Functional<double>
{
public:
  SH_smooth_energy_tet(const mati_t &tets, const matd_t &nods, const double w)
      : tets_(tets), w_(w), dim_(3*tets.size(2)) {
    matd_t volume = zeros<double>(tets.size(2), 1); {
      #pragma omp parallel for
      for (size_t i = 0; i < tets.size(2); ++i) {
        matd_t Ds = nods(colon(), tets(colon(1, 3), i))-nods(colon(), tets(0, i))*ones<double>(1, 3);
        volume[i] = fabs(det(Ds))/6.0;
      }
    }
    
    shared_ptr<face2tet_adjacent> f2t(face2tet_adjacent::create(tets));
    vector<size_t> buffer;
    for (size_t i = 0; i < f2t->face2tet_.size(); ++i) {
      const pair<size_t, size_t> ts = f2t->face2tet_[i];
      if ( !f2t->is_outside_face(ts) ) {
        buffer.push_back(ts.first);
        buffer.push_back(ts.second);
      }
    }
    adjt_.resize(2, buffer.size()/2);
    std::copy(buffer.begin(), buffer.end(), adjt_.begin());

    stiff_.resize(adjt_.size(2)); {
      #pragma omp parallel for
      for (size_t i = 0; i < stiff_.size(); ++i) {
        const double dist = norm(
            nods(colon(), tets(colon(), adjt_(0, i)))*ones<double>(4, 1)/4.0-
            nods(colon(), tets(colon(), adjt_(1, i)))*ones<double>(4, 1)/4.0);
        stiff_[i] = (volume[adjt_(0, i)]+volume[adjt_(1, i)])/(dist*dist);
      }
      double total = sum(stiff_);
      stiff_ /= total;
    }
  }
  size_t Nx() const {
    return dim_;
  }
  int Val(const double *abc, double *val) const {
    itr_matrix<const double *> ABC(3, dim_/3, abc);
    matd_t abcs = zeros<double>(3, 2); double value = 0;
    for (size_t i = 0; i < adjt_.size(2); ++i) {
      abcs = ABC(colon(), adjt_(colon(), i));
      cubic_sym_smooth_tet_(&value, &abcs[0], &stiff_[i]);
      *val += w_*value;
    }
    return 0;
  }
  int Gra(const double *abc, double *gra) const {
    itr_matrix<const double *> ABC(3, dim_/3, abc);
    itr_matrix<double *> G(3, dim_/3, gra);
    matd_t abcs = zeros<double>(3, 2), g = zeros<double>(3, 2);
    for (size_t i = 0; i < adjt_.size(2); ++i) {
      abcs = ABC(colon(), adjt_(colon(), i));
      cubic_sym_smooth_tet_jac_(&g[0], &abcs[0], &stiff_[i]);
      G(colon(), adjt_(colon(), i)) += w_*g;
    }
    return 0;
  }
  int Hes(const double *abc, vector<Triplet<double>> *hes) const {
    return __LINE__;
  }
  int ValSH(const double *f, double *val) const {
    itr_matrix<const double *> Fs(9, dim_/3, f);
    matd_t fs = zeros<double>(9, 2), dif = zeros<double>(9, 1);
    for (size_t i = 0; i < adjt_.size(2); ++i) {
      fs = Fs(colon(), adjt_(colon(), i));
      dif = fs(colon(), 0)-fs(colon(), 1);
      *val += w_*stiff_[i]*dot(dif, dif);
    }
    return 0;
  }
  int GraSH(const double *f, double *gra) const {
    itr_matrix<const double *> Fs(9, dim_/3, f);
    itr_matrix<double *> G(9, dim_/3, gra);
    matd_t fs = zeros<double>(9, 2);
    for (size_t i = 0; i < adjt_.size(2); ++i) {
      fs = Fs(colon(), adjt_(colon(), i));
      G(colon(), adjt_(0, i)) += 2.0*w_*stiff_[i]*(fs(colon(), 0)-fs(colon(), 1));
      G(colon(), adjt_(1, i)) += 2.0*w_*stiff_[i]*(fs(colon(), 1)-fs(colon(), 0));
    }
    return 0;
  }
  int HesSH(const double *f, vector<Triplet<double>> *hes) const {
    for (size_t i = 0; i < adjt_.size(2); ++i) {
      const size_t l = adjt_(0, i), r = adjt_(1, i);
      const double entry = 2.0*w_*stiff_[i];
      add_diag_block<double, 9>(l, l, entry, hes);
      add_diag_block<double, 9>(l, r, -entry, hes);
      add_diag_block<double, 9>(r, l, -entry, hes);
      add_diag_block<double, 9>(r, r, entry, hes);
    }
    return 0;
  }
private:
  const mati_t &tets_;
  const double w_;
  const size_t dim_;

  mati_t adjt_;
  matd_t stiff_;
};

class SH_align_energy_tet : public Functional<double>
{
public:
  SH_align_energy_tet(const mati_t &tets, const matd_t &nods, const double w)
      : tets_(tets), w_(w), dim_(3*tets.size(2)) {
    shared_ptr<face2tet_adjacent> f2t(face2tet_adjacent::create(tets));
    mati_t surf; bool check_order = true;
    jtf::mesh::get_outside_face(*f2t, surf, check_order, &nods);
  
    stiff_.resize(surf.size(2)); {
      #pragma omp parallel for
      for (size_t i = 0; i < surf.size(2); ++i)
        stiff_[i] = jtf::mesh::cal_face_area(nods(colon(), surf(colon(), i)));
      double sum_area = sum(stiff_);
      stiff_ /= sum_area;
    }
    
    adjt_.resize(surf.size(2));
    zyz_.resize(3, surf.size(2)); {
      #pragma omp parallel for
      for (size_t i = 0; i < surf.size(2); ++i) {
        matd_t n = zeros<double>(3, 1);
        jtf::mesh::cal_face_normal(nods(colon(), surf(colon(), i)), n);
        normal2zyz(&n[0], &zyz_(0, i));

        const pair<size_t, size_t> t = f2t->query(surf(0, i), surf(1, i), surf(2, i));
        adjt_[i] = (t.first == -1) ? t.second : t.first;
      }
    }
  }
  size_t Nx() const {
    return dim_;
  }
  int Val(const double *abc, double *val) const {
    itr_matrix<const double *> ABC(3, dim_/3, abc);
    double value = 0;
    for (size_t i = 0; i < adjt_.size(); ++i) {
      const size_t tid = adjt_[i];
      cubic_sym_align_(&value, &ABC(0, tid), &zyz_(0, i), &stiff_[i]);
      *val += w_*value;
    }
    return 0;
  }
  int Gra(const double *abc, double *gra) const {
    itr_matrix<const double *> ABC(3, dim_/3, abc);
    itr_matrix<double *> G(3, dim_/3, gra);
    matd_t g = zeros<double>(3, 1);
    for (size_t i = 0; i < adjt_.size(); ++i) {
      const size_t tid = adjt_[i];
      cubic_sym_align_jac_(&g[0], &ABC(0, tid), &zyz_(0, i), &stiff_[i]);
      G(colon(), tid) += w_*g;
    }
    return 0;
  }
  int Hes(const double *abc, vector<Triplet<double>> *hes) const {
    return __LINE__;
  }
  int ValSH(const double *f, double *val) const {
    itr_matrix<const double *> Fs(9, dim_/3, f);
    double value = 0;
    for (size_t i = 0; i < adjt_.size(); ++i) {
      const size_t tid = adjt_[i];
      cubic_align_sh_coef_(&value, &Fs(0, tid), &zyz_(0, i), &stiff_[i]);
      *val += w_*value;
    }
    return 0;
  }
  int GraSH(const double *f, double *gra) const {
    itr_matrix<const double *> Fs(9, dim_/3, f);
    itr_matrix<double *> G(9, dim_/3, gra);
    matd_t g = zeros<double>(9, 1);
    for (size_t i = 0; i < adjt_.size(); ++i) {
      const size_t tid = adjt_[i];
      cubic_align_sh_coef_jac_(&g[0], &Fs(0, tid), &zyz_(0, i), &stiff_[i]);
      G(colon(), tid) += w_*g;
    }
    return 0;
  }
  int HesSH(const double *f, vector<Triplet<double>> *hes) const {
    matd_t H = zeros<double>(9, 9);
    for (size_t i = 0; i < adjt_.size(); ++i) {
      const size_t tid = adjt_[i];
      cubic_align_sh_coef_hes_(&H[0], nullptr, &zyz_(0, i), &stiff_[i]);
      for (size_t p = 0; p < 9; ++p) {
        for (size_t q = 0; q < 9; ++q) {
          const size_t I = 9*tid+p;
          const size_t J = 9*tid+q;
          hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
        }
      }
    }
    return 0;
  }
private:
  const mati_t &tets_;
  const double w_;
  const size_t dim_;

  mati_t adjt_;
  matd_t stiff_;
  matd_t zyz_;
};

//===============================================================================

SH_smooth_energy::SH_smooth_energy(const mati_t &tets, const matd_t &nods, const double w)
    : tets_(tets), w_(w), dim_(3*nods.size(2)) {
  vol_ = zeros<double>(1, tets_.size(2)); {
    #pragma omp parallel for
    for (size_t i = 0; i < tets_.size(2); ++i) {
      matd_t Ds = nods(colon(), tets_(colon(1, 3), i))-nods(colon(), tets_(0, i))*ones<double>(1, 3);
      vol_[i] = fabs(det(Ds))/6.0;
    }
  }

  w_ /= std::cbrt(sum(vol_));
  
  basis_ = zeros<double>(12, tets_.size(2)); {
    #pragma omp parallel for
    for (size_t i = 0; i < tets_.size(2); ++i) {
      matd_t vert = nods(colon(), tets_(colon(), i));
      calc_tet_linear_basis_grad(&vert[0], &basis_(0, i));
    }
  }
}

size_t SH_smooth_energy::Nx() const {
  return dim_;
}

int SH_smooth_energy::Val(const double *abc, double *val) const {
  itr_matrix<const double *> ABC(3, dim_/3, abc);
  matd_t abcs = zeros<double>(3, 4); double value = 0;
  for (size_t i = 0; i < tets_.size(2); ++i) {
    abcs = ABC(colon(), tets_(colon(), i));
    cubic_sym_smooth_(&value, &abcs[0], &basis_(0, i), &vol_[i]);
    *val += w_*value;
  }
  return 0;
}

int SH_smooth_energy::Gra(const double *abc, double *gra) const {
  itr_matrix<const double *> ABC(3, dim_/3, abc);
  itr_matrix<double *> G(3, dim_/3, gra);
  matd_t abcs = zeros<double>(3, 4), g = zeros<double>(3, 4);
  for (size_t i = 0; i < tets_.size(2); ++i) {
    abcs = ABC(colon(), tets_(colon(), i));
    cubic_sym_smooth_jac_(&g[0], &abcs[0], &basis_(0, i), &vol_[i]);
    G(colon(), tets_(colon(), i)) += w_*g;
  }
  return 0;
}

int SH_smooth_energy::ValSH(const double *f, double *val) const {
  itr_matrix<const double *> F(9, dim_/3, f);
  matd_t fs = zeros<double>(9, 4); double value = 0;
  for (size_t i = 0; i < tets_.size(2); ++i) {
    fs = F(colon(), tets_(colon(), i));
    cubic_smooth_sh_coef_(&value, &fs[0], &basis_(0, i), &vol_[i]);
    *val += w_*value;
  }
  return 0;
}

int SH_smooth_energy::GraSH(const double *f, double *gra) const {
  itr_matrix<const double *> F(9, dim_/3, f);
  itr_matrix<double *> G(9, dim_/3, gra);
  matd_t fs = zeros<double>(9, 4), g = zeros<double>(9, 4);
  for (size_t i = 0; i < tets_.size(2); ++i) {
    fs = F(colon(), tets_(colon(), i));
    cubic_smooth_sh_coef_jac_(&g[0], &fs[0], &basis_(0, i), &vol_[i]);
    G(colon(), tets_(colon(), i)) += w_*g;
  }
  return 0;
}

int SH_smooth_energy::HesSH(const double *f, vector<Triplet<double>> *hes) const {
  matd_t H = zeros<double>(36, 36);
  for (size_t i = 0; i < tets_.size(2); ++i) {
    cubic_smooth_sh_coef_hes_(&H[0], nullptr, &basis_(0, i), &vol_[i]);
    for (size_t p = 0; p < 36; ++p) {
      for (size_t q = 0; q < 36; ++q) {
        if ( H(p, q) != 0.0 ) {
          const size_t I = 9*tets_(p/9, i)+p%9;
          const size_t J = 9*tets_(q/9, i)+q%9;
          hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
        }
      }
    }
  }
  return 0;
}
//===============================================================================
/// @return zyz angle to align normal and axis z
static inline void normal2zyz(const double *n, double *zyz) {
  zyz[0] = -atan2(n[1], n[0]);
  zyz[1] = -acos(n[2]);
  zyz[2] = 0;
}

SH_align_energy::SH_align_energy(const mati_t &tets, const matd_t &nods, const double w)
    : tets_(tets), w_(w), dim_(3*nods.size(2)) {
  shared_ptr<face2tet_adjacent> f2t(face2tet_adjacent::create(tets));
  bool check_order = true;
  jtf::mesh::get_outside_face(*f2t, surf_, check_order, &nods);
  cout << "[INFO] surface elements: " << surf_.size(2) << endl;
  
  area_.resize(1, surf_.size(2)); {
    #pragma omp parallel for
    for (size_t i = 0; i < surf_.size(2); ++i)
      area_[i] = jtf::mesh::cal_face_area(nods(colon(), surf_(colon(), i)));
  }

  w_ /= sum(area_);

  zyz_.resize(3, surf_.size(2)); {
    #pragma omp parallel for
    for (size_t i = 0; i < surf_.size(2); ++i) {
      matd_t n = zeros<double>(3, 1);
      jtf::mesh::cal_face_normal(nods(colon(), surf_(colon(), i)), n);
      normal2zyz(&n[0], &zyz_(0, i));
    }
  }
}

size_t SH_align_energy::Nx() const {
  return dim_;
}

int SH_align_energy::Val(const double *abc, double *val) const {
  itr_matrix<const double *> ABC(3, dim_/3, abc);
  double value = 0;
  for (size_t i = 0; i < surf_.size(2); ++i) {
    for (size_t j = 0; j < surf_.size(1); ++j) {
      const size_t idx = surf_(j, i);
      cubic_sym_align_(&value, &ABC(0, idx), &zyz_(0, i), &area_[i]);
      *val += w_*value;
    }
  }
  return 0;
}

int SH_align_energy::Gra(const double *abc, double *gra) const {
  itr_matrix<const double *> ABC(3, dim_/3, abc);
  itr_matrix<double *> G(3, dim_/3, gra);
  matd_t g = zeros<double>(3, 1);
  for (size_t i = 0; i < surf_.size(2); ++i) {
    for (size_t j = 0; j < surf_.size(1); ++j) {
      const double idx = surf_(j, i);
      cubic_sym_align_jac_(&g[0], &ABC(0, idx), &zyz_(0, i), &area_[i]);
      G(colon(), idx) += w_*g;
    }
  }
  return 0;
}

int SH_align_energy::ValSH(const double *f, double *val) const {
  itr_matrix<const double *> Fs(9, dim_/3, f);
  double value = 0;
  for (size_t i = 0; i < surf_.size(2); ++i) {
    for (size_t j = 0; j < surf_.size(1); ++j) {
      const size_t idx = surf_(j, i);
      cubic_align_sh_coef_(&value, &Fs(0, idx), &zyz_(0, i), &area_[i]);
      *val += w_*value;
    }
  }
  return 0;
}

int SH_align_energy::GraSH(const double *f, double *gra) const {
  itr_matrix<const double *> Fs(9, dim_/3, f);
  itr_matrix<double *> G(9, dim_/3, gra);
  matd_t g = zeros<double>(9, 1);
  for (size_t i = 0; i < surf_.size(2); ++i) {
    for (size_t j = 0; j < surf_.size(1); ++j) {
      const size_t idx = surf_(j, i);
      cubic_align_sh_coef_jac_(&g[0], &Fs(0, idx), &zyz_(0, i), &area_[i]);
      G(colon(), idx) += w_*g;
    }
  }
  return 0;
}

int SH_align_energy::HesSH(const double *f, vector<Triplet<double>> *hes) const {
  matd_t H = zeros<double>(9, 9);
  for (size_t i = 0; i < surf_.size(2); ++i) {
    for (size_t j = 0; j < surf_.size(1); ++j) {
      const size_t idx = surf_(j, i);
      cubic_align_sh_coef_hes_(&H[0], nullptr, &zyz_(0, i), &area_[i]);
      for (size_t p = 0; p < 9; ++p) {
        for (size_t q = 0; q < 9; ++q) {
          if ( H(p, q) != 0.0 ) {
            const size_t I = 9*idx+p;
            const size_t J = 9*idx+q;
            hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
          }
        }
      }
    }
  }
  return 0;
}
//===============================================================================
int cross_frame_opt::init(const mati_t &tets, const matd_t &nods, const cross_frame_args &args) {
  args_ = args;

  int success = 0;
  buffer_.push_back(make_shared<SH_smooth_energy_tet>(tets, nods, args_.ws));
  buffer_.push_back(make_shared<SH_align_energy_tet>(tets, nods, args_.wa));
  try {
    energy_ = make_shared<energy_t<double>>(buffer_);
  } catch ( exception &e ) {
    cerr << "[Exception] " << e.what() << endl;
    success = 1;
  }
  return success;
}

cross_frame_opt* cross_frame_opt::create(const mati_t &tets, const matd_t &nods, const cross_frame_args &args) {
  cross_frame_opt *handle = new cross_frame_opt;
  if ( handle->init(tets, nods, args) )
    return nullptr;
  return handle;
}

int cross_frame_opt::solve_laplacian(VectorXd &Fs) const {
  const size_t dim = 3*energy_->Nx();
  Fs = VectorXd::Zero(dim);
  cout << "\t@linear solve dimension: " << dim << endl << endl;

  auto fs = dynamic_pointer_cast<SH_smooth_energy_tet>(buffer_[0]);
  auto fa = dynamic_pointer_cast<SH_align_energy_tet>(buffer_[1]);

  double prev_vs = 0, prev_va = 0; {
    fs->ValSH(Fs.data(), &prev_vs);
    fa->ValSH(Fs.data(), &prev_va);
    cout << "\t@prev smoothness energy: " << prev_vs << endl;
    cout << "\t@prev alignment energy: " << prev_va << endl << endl;
  }

  VectorXd g = VectorXd::Zero(dim); {
    fs->GraSH(Fs.data(), g.data());
    fa->GraSH(Fs.data(), g.data());
    g *= -1;
    cout << "\t@prev grad norm: " << g.norm() << endl << endl;
  }
  
  SparseMatrix<double> H(dim, dim); {
    vector<Triplet<double>> trips;
    fs->HesSH(nullptr, &trips);
    fa->HesSH(nullptr, &trips);
    H.setFromTriplets(trips.begin(), trips.end());
    H.makeCompressed();
  }

#ifdef USE_CHOLMOD
  CholmodSimplicialLLT<SparseMatrix<double>> solver;
  solver.compute(H);
  ASSERT(solver.info() == Success);
  VectorXd dx = solver.solve(g);
  ASSERT(solver.info() == Success);
#else
  boost::property_tree::ptree pt;
  pt.put("linear_solver/type.value", "PETsc");
  shared_ptr<linear_solver> solver
      (linear_solver::create(H.valuePtr(), H.innerIndexPtr(), H.outerIndexPtr(), H.nonZeros(), dim, dim, pt));
  VectorXd dx = VectorXd::Zero(dim);
  solver->solve(g.data(), dx.data(), dim, pt);
#endif
  
  Fs += dx;
  
  double post_vs = 0, post_va = 0; {
    fs->ValSH(Fs.data(), &post_vs);
    fa->ValSH(Fs.data(), &post_va);
    cout << "\t@post smoothness energy: " << post_vs << endl;
    cout << "\t@post alignment energy: " << post_va << endl << endl;
  }

  VectorXd gs = VectorXd::Zero(dim); {
    fs->GraSH(Fs.data(), gs.data());
    fa->GraSH(Fs.data(), gs.data());
    cout << "\t@post grad norm: " << gs.norm() << endl << endl;
  }

  cout << "\t@SOLUTION NORM: " << Fs.norm() << endl;
  return 0;
}

int cross_frame_opt::solve_initial_frames(const VectorXd &Fs, VectorXd &abc) const {
  abc = VectorXd::Zero(energy_->Nx());
  
  #pragma omp parallel for
  for (size_t i = 0; i < abc.size()/3; ++i) {
    // double prev_res = 0; 
    // sh_residual_(&prev_res, &abc[3*i], &Fs[9*i]);
    
    sh_to_zyz(&Fs[9*i], &abc[3*i], 1000);

    // double post_res = 0;
    // sh_residual_(&post_res, &abc[3*i], &Fs[9*i]);
  }

  return 0;
}

int cross_frame_opt::optimize_frames(VectorXd &abc) const {
  const double epsf = args_.epsf, epsx = 0;
  const size_t maxits = args_.maxits;

  {
    double vs = 0, va = 0;
    buffer_[0]->Val(abc.data(), &vs);
    buffer_[1]->Val(abc.data(), &va);
    cout << "\t@prev smoothness energy: " << vs << endl;
    cout << "\t@prev alignment energy: " << va << endl;
    VectorXd g = VectorXd::Zero(energy_->Nx());
    energy_->Gra(abc.data(), g.data());
    cout << "\t@prev gradient norm: " << g.norm() << endl << endl;
  }
  
  lbfgs_solve(energy_, abc.data(), abc.size(), epsf, epsx, maxits);

  {
    double vs = 0, va = 0;
    buffer_[0]->Val(abc.data(), &vs);
    buffer_[1]->Val(abc.data(), &va);
    cout << "\t@post smoothness energy: " << vs << endl;
    cout << "\t@post alignment energy: " << va << endl;
    VectorXd g = VectorXd::Zero(energy_->Nx());
    energy_->Gra(abc.data(), g.data());
    cout << "\t@post gradient norm: " << g.norm() << endl << endl;
  }

  return 0;
}

}
