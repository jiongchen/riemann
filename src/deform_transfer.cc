#include "deform_transfer.h"

#include <jtflib/mesh/io.h>

#include "vtk.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;
using namespace surfparam;

namespace geom_deform {

extern "C" {

void deform_energy_unit_();
void deform_energy_unit_jac_();
void deform_energy_unit_hes_();

void unit_smooth_energy_();
void unit_smooth_energy_jac_();
void unit_smooth_energy_hes_();

}

class deform_energy : public Functional<double>
{

};

class smooth_energy : public Functional<double>
{

};

class identity_energy : public Functional<double>
{

};

class closest_energy : public Functional<double>
{

};

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

}

int deform_transfer::save_deformed_target_mesh(const char *filename) const {

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

}
