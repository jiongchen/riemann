#include "volume_frame.h"

#include <zjucad/matrix/itr_matrix.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>
#include <Eigen/Eigen>
#include <Eigen/CholmodSupport>

#include "config.h"
#include "grad_operator.h"

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
  void cubic_smooth_sh_coef_hes_(double *hes, const double *F, const double *CR, const double *vol);;

  void cubic_align_sh_coef_(double *val, const double *F, const double *Rnz, const double *area);
  void cubic_align_sh_coef_jac_(double *jac, const double *F, const double *Rnz, const double *area);
  void cubic_align_sh_coef_hes_(double *hes, const double *F, const double *Rnz, const double *area);
                                
  void sh_residual_(double *val, const double *abc, const double *f);
  void sh_residual_jac_(double *jac, const double *abc, const double *f);
  void sh_residual_hes_(double *hes, const double *abc, const double *f);

}

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
  matd_t abcs = zeros<double>(3, 4);
  for (size_t i = 0; i < tets_.size(2); ++i) {
    abcs = ABC(colon(), tets_(colon(), i));
    double value = 0;
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
  matd_t fs = zeros<double>(9, 4);
  double value = 0;
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
        if ( fabs(H(p, q)) >= 1e-12 ) {
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
  matd_t abcs = zeros<double>(3, 1);
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
  matd_t abcs = zeros<double>(3, 1), g = zeros<double>(3, 1);
  for (size_t i = 0; i < surf_.size(2); ++i) {
    for (size_t j = 0; j < 3; ++j) {
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
          const size_t I = 9*idx+p;
          const size_t J = 9*idx+q;
          hes->push_back(Triplet<double>(I, J, w_*H(p, q)));
        }
      }
    }
  }
  return 0;
}
//===============================================================================
// static void reconstruct_zyz(double *abc, const double *f, const size_t maxiter) {
//   matd_t g = zeros<double>(3, 1);
//   matd_t H = zeros<double>(3, 3);
//   for (size_t i = 0; i < maxiter; ++i) {
//     sh_residual_jac_(&g[0], abc, f);
//     sh_residual_hes_(&H[0], abc, f);
//     if ( inv(H) ) {
//       cerr << "\t@fail to inverse" << endl;
//       return;
//     }
//     itr_matrix<double *>(3, 1, abc) += -H*g;
//   }
// }

int cross_frame_opt::init(const mati_t &tets, const matd_t &nods) {
  int success = 0;
  buffer_.push_back(make_shared<SH_smooth_energy>(tets, nods, 1e0));
  buffer_.push_back(make_shared<SH_align_energy>(tets, nods, 1e3));
  try {
    energy_ = make_shared<energy_t<double>>(buffer_);
  } catch ( exception &e ) {
    cerr << "[Exception] " << e.what() << endl;
    success = 1;
  }
  return success;
}

cross_frame_opt* cross_frame_opt::create(const mati_t &tets, const matd_t &nods) {
  cross_frame_opt *handle = new cross_frame_opt;
  if ( handle->init(tets, nods) )
    return nullptr;
  return handle;
}

int cross_frame_opt::solve_smooth_sh_coeffs(VectorXd &Fs) const {
  const size_t dim = Fs.size();
  auto f1 = dynamic_pointer_cast<SH_smooth_energy>(buffer_[0]);
  auto f2 = dynamic_pointer_cast<SH_align_energy>(buffer_[1]);
  double prev_value = 0; {
    f1->ValSH(Fs.data(), &prev_value);
    f2->ValSH(Fs.data(), &prev_value);
    cout << "[Info] prev energy value: " << prev_value << endl;
  }
  VectorXd g = VectorXd::Zero(dim); {
    f1->GraSH(Fs.data(), g.data());
    f2->GraSH(Fs.data(), g.data());
  }
  SparseMatrix<double> H(dim, dim); {
    vector<Triplet<double>> trips;
    f1->HesSH(Fs.data(), &trips);
    f2->HesSH(Fs.data(), &trips);
    H.reserve(trips.size());
    H.setFromTriplets(trips.begin(), trips.end());
  }
  CholmodSimplicialLDLT<SparseMatrix<double>> solver;
  solver.compute(H);
  ASSERT(solver.info() == Success);
  VectorXd dx = -solver.solve(g);
  ASSERT(solver.info() == Success);
  Fs += dx;
  double post_value = 0; {
    f1->ValSH(Fs.data(), &post_value);
    f2->ValSH(Fs.data(), &post_value);
    cout << "[Info] post value: " << post_value << endl;
  }
  return 0;
}

int cross_frame_opt::solve_initial_frames(VectorXd &abc) const {
  return 0;
}

int cross_frame_opt::optimize_frames(VectorXd &abc) const {
  return 0;
}

}
