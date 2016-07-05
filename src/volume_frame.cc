#include "volume_frame.h"

#include <zjucad/matrix/itr_matrix.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/util.h>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>

using namespace std;
using namespace zjucad::matrix;
using namespace jtf::mesh;

namespace riemann {

extern "C" {

void cubic_sym_smooth_(double *val, const double *abc, const double *CR, const double *vol);
void cubic_sym_smooth_jac_(double *jac, const double *abc, const double *CR, const double *vol);

void cubic_sym_align_(double *val, const double *abc, const double *rnz, const double *area);
void cubic_sym_align_jac_(double *jac, const double *abc, const double *rnz, const double *area);

void sh_residual_(double *val, const double *abc, const double *f);
void sh_residual_jac_(double *jac, const double *abc, const double *f);

}

//===============================================================================
/// @return Crouzeix-Raviart basis
static void tet_cr_basis(const double *x, double *basis) {
  
}

SH_smooth_energy::SH_smooth_energy(const mati_t &tets, const matd_t &nods, const double w)
    : tets_(tets), nods_(nods), w_(w) {
  shared_ptr<face2tet_adjacent> f2t(face2tet_adjacent::create(tets));
  dim_ = 3*f2t->faces_.size();

  vol_.resize(1, tets_.size(2)); {
    for (size_t i = 0; i < tets_.size(2); ++i) {
      matd_t Ds = nods_(colon(), tets_(colon(1, 3), i))-nods_(colon(), tets_(0, i))*ones<double>(1, 3);
      vol_[i] = fabs(det(Ds))/6.0;
    }
  }

  w_ /= sum(vol_);
  
  CR_.resize(12, tets_.size(2)); {
    // TODO
  }
  
  tet_face_.resize(4, tets_.size(2)); {
    for (size_t i = 0; i < tets_.size(2); ++i) 
      for (size_t j = 0; j < 4; ++j)
        tet_face_(j, i) = f2t->get_face_idx(tets_(j, i), tets_((j+1)%4, i), tets_((j+2)%4, i));
  }
}

size_t SH_smooth_energy::Nx() const {
  return dim_;
}

int SH_smooth_energy::Val(const double *abc, double *val) const {
  itr_matrix<const double *> ABC(3, dim_/3, abc);
  matd_t abcs(3, 4);
  for (size_t i = 0; i < tet_face_.size(2); ++i) {
    abcs = ABC(colon(), tet_face_(colon(), i));
    double value = 0;
    cubic_sym_smooth_(&value, &abcs[0], &CR_(0, i), &vol_[i]);
    *val += w_*value;
  }
  return 0;
}

int SH_smooth_energy::Gra(const double *abc, double *gra) const {
  itr_matrix<const double *> ABC(3, dim_/3, abc);
  itr_matrix<double *> G(3, dim_/3, gra);
  matd_t abcs(3, 4), g(3, 4);
  for (size_t i = 0; i < tet_face_.size(2); ++i) {
    abcs = ABC(colon(), tet_face_(colon(), i));
    cubic_sym_smooth_jac_(&g[0], &abcs[0], &CR_(0, i), &vol_[i]);
    G(colon(), tet_face_(colon(), i)) += w_*g;
  }
  return 0;
}
//===============================================================================
/// @return zyz angle to align normal and axis z
static void normal2zyz(const double *n, double *zyz) {
  
}

SH_align_energy::SH_align_energy(const mati_t &tets, const matd_t &nods, const double w)
    : tets_(tets), nods_(nods), w_(w) {
  shared_ptr<face2tet_adjacent> f2t(face2tet_adjacent::create(tets));
  dim_ = 3*f2t->faces_.size();
  
  mati_t surf; bool check_order = true;
  jtf::mesh::get_outside_face(*f2t, surf, check_order, &nods);

  area_.resize(1, surf.size(2)); {
    for (size_t i = 0; i < surf.size(2); ++i)
      area_[i] = jtf::mesh::cal_face_area(nods_(surf(colon(), i)));
  }

  w_ /= sum(area_);

  zyz_.resize(3, surf.size(2)); {
    matd_t n(3, 1);
    for (size_t i = 0; i < surf.size(2); ++i) {
      jtf::mesh::cal_face_normal(nods_(colon(), surf(colon(), i)), n);
      normal2zyz(&n[0], &zyz_(0, i));
    }
  }

  bnd_face_.resize(1, surf.size(2)); {
    for (size_t i = 0; i < surf.size(2); ++i)
      bnd_face_[i] = f2t->get_face_idx(surf(0, i), surf(1, i), surf(2, i));
  }
}

size_t SH_align_energy::Nx() const {
  return dim_;
}

int SH_align_energy::Val(const double *abc, double *val) const {
  itr_matrix<const double *> ABC(3, dim_/3, abc);
  double value = 0;
  for (size_t i = 0; i < bnd_face_.size(); ++i) {
    const double idx = bnd_face_[i];
    cubic_sym_align_(&value, &ABC(0, idx), &zyz_(0, i), &area_[i]);
    *val += w_*value;
  }
  return 0;
}

int SH_align_energy::Gra(const double *abc, double *gra) const {
  itr_matrix<const double *> ABC(3, dim_/3, abc);
  itr_matrix<double *> G(3, dim_/3, gra);
  matd_t jac = zeros<double>(3, 1);
  for (size_t i = 0; i < bnd_face_.size(); ++i) {
    const double idx = bnd_face_[i];
    cubic_sym_align_jac_(&jac[0], &ABC(0, idx), &zyz_(0, i), &area_[i]);
    G(colon(), idx) += w_*jac;
  }
  return 0;
}
//===============================================================================
//===============================================================================
}
