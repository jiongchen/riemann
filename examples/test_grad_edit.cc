#include <iostream>
#include <jtflib/mesh/io.h>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <Eigen/Dense>

#include "src/gradient_deform.h"
#include "src/write_vtk.h"

using namespace std;
using namespace riemann;
using namespace Eigen;

static int read_fixed_vert(const char *filename, vector<size_t> &vert) {
  ifstream ifs(filename);
  if ( ifs.fail() ) {
    cerr << "[INFO] cannot open " << filename << endl;
    return __LINE__;
  }
  size_t id;
  while ( ifs >> id ) {
    vert.push_back(id);
  }
  cout << "[INFO] constrained verts: " << vert.size() << endl;
  ifs.close();
  return 0;
}

namespace test_grad_edit {
struct argument {
  string input_obj_file;
  string fix_vert_file;
  string edit_vert_file;
  string output_folder;
};
}

int main(int argc, char *argv[])
{
  if ( argc != 4 ) {
    cerr << "# Usage: " << argv[0] << " model.obj fix.fv edit.fv\n";
    return __LINE__;
  }
  boost::filesystem::create_directory("./grad_edit");

  mati_t tris; matd_t nods;
  jtf::mesh::load_obj(argv[1], tris, nods);
  vector<size_t> fix_vert, edit_vert;
  read_fixed_vert(argv[2], fix_vert);
  read_fixed_vert(argv[3], edit_vert);

  gradient_field_deform dfm(tris, nods);
  dfm.set_fixed_verts(fix_vert);
  dfm.set_edited_verts(edit_vert);

  Quaterniond Q; {
    Vector3d axis(0, -1, 0);
    double angle = 0.75*M_PI;
    Q = Quaterniond(AngleAxisd(angle, axis));
  }
  dfm.prescribe_uniform_transform(Q, 1.0);
//  dfm.precompute();

  mati_t rtris; matd_t rnods;
  {
    dfm.solve_harmonic_field();
    dfm.propogate_transform();
    dfm.interp_face_transform();
    dfm.rotate_surf_piecewise(rtris, rnods);
    dfm.update_gradient_field();
    dfm.prefactorize();
  }

  dfm.deform(&nods[0]);

  MatrixXd dG = dfm.Gxyz_;
  dG *= 0.1;

  char outfile[256];
  sprintf(outfile, "./grad_edit/harmonic.vtk");
  draw_vert_value_to_vtk(outfile, &nods[0], nods.size(2), &tris[0], tris.size(2), &dfm.hf_[0]);
  sprintf(outfile, "./grad_edit/gradX.vtk");
  draw_face_direct_field(outfile, &nods[0], nods.size(2), &tris[0], tris.size(2), &dG(0, 0));
  sprintf(outfile, "./grad_edit/gradY.vtk");
  draw_face_direct_field(outfile, &nods[0], nods.size(2), &tris[0], tris.size(2), &dG(0, 1));
  sprintf(outfile, "./grad_edit/gradZ.vtk");
  draw_face_direct_field(outfile, &nods[0], nods.size(2), &tris[0], tris.size(2), &dG(0, 2));
  sprintf(outfile, "./grad_edit/rotated_fragment.obj");
  jtf::mesh::save_obj(outfile, rtris, rnods);
  sprintf(outfile, "./grad_edit/deformed.obj");
  jtf::mesh::save_obj(outfile, tris, nods);

  cout << "[INFO] done\n";
  return 0;
}
