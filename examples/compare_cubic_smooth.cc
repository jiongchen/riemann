#include <iostream>
#include <fstream>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <zjucad/ptree/ptree.h>

#include "src/def.h"
#include "src/volume_frame.h"
#include "src/vtk.h"
#include "src/write_vtk.h"

using namespace std;
using namespace riemann;
using namespace Eigen;
using namespace zjucad::matrix;

static int write_tet_zyz(const char *filename, const double *zyz, const size_t elem_num) {
  ofstream ofs(filename);
  if ( ofs.fail() ) {
    cerr << "[Error] can not write to " << filename << endl;
    return __LINE__;
  }
  ofs << "3 " << elem_num << endl;
  for (size_t i = 0; i < elem_num; ++i)
    ofs << zyz[3*i+0] << " " << zyz[3*i+1] << " " << zyz[3*i+2] << endl;
  return 0;
}

static int write_tet_ff(const char *filename, const double *ff, const size_t elem_num) {
  ofstream ofs(filename);
  if ( ofs.fail() ) {
    cerr << "[Error] can not write to " << filename << endl;
    return __LINE__;
  }
  ofs << "9 " << elem_num << endl;
  for (size_t i = 0; i < elem_num; ++i)
    ofs << ff[9*i+0] << " " << ff[9*i+1] << " " << ff[9*i+2] << " "
        << ff[9*i+3] << " " << ff[9*i+4] << " " << ff[9*i+5] << " "
        << ff[9*i+6] << " " << ff[9*i+7] << " " << ff[9*i+8] << endl;
  return 0;
}

int main(int argc, char *argv[])
{
  boost::property_tree::ptree pt;
  zjucad::read_cmdline(argc, argv, pt);

  string out_folder = pt.get<string>("out_dir.value");
  
  mati_t tets; matd_t nods;
  jtf::mesh::tet_mesh_read_from_vtk(pt.get<string>("mesh.value").c_str(), &nods, &tets);
  {    
    string outfile = out_folder+string("/tet.vtk");
    ofstream ofs(outfile);
    tet2vtk(ofs, &nods[0], nods.size(2), &tets[0], tets.size(2));
    ofs.close();
  }

  // GET INITAL VALUE
  boost::property_tree::ptree ptt; {
    ptt.put("weight.smooth.value", 1e-1);
    ptt.put("weight.align.value", 1e3);
    ptt.put("lbfgs.epsf.value", 1e-8);
    ptt.put("lbfgs.maxits.value", 2);
  }

  shared_ptr<cross_frame_opt> frame_opt = make_shared<cross_frame_opt>(tets, nods, ptt);

  cout << "[INFO] solve Laplacian" << endl;
  VectorXd Fs;
  frame_opt->solve_laplacian(Fs);

  cout << "[INFO] solve for initial zyz angles" << endl;
  VectorXd abc;
  frame_opt->solve_initial_frames(Fs, abc);

  cout << "[INFO] optimize frames" << endl;
  frame_opt->optimize_frames(abc);

  // write initial zyz
  string init_zyz_file = out_folder+string("/init_zyz.txt");
  write_tet_zyz(init_zyz_file.c_str(), abc.data(), abc.size()/3);
  
  // OPTIMZE SMOOTHNESS
  cout << "[INFO] init frame smoother" << endl;
  shared_ptr<frame_smoother> smoother = make_shared<frame_smoother>(tets, nods, pt);

  if ( pt.get<string>("sm_type.value") == "SH" ) {
    
    cout << "[INFO] smooth internal frames via SH" << endl;
    smoother->smoothSH(abc);

    string zyz_file = out_folder+string("/frames.txt");
    write_tet_zyz(zyz_file.c_str(), abc.data(), abc.size()/3);

  } else if ( pt.get<string>("sm_type.value") == "L1" ) {

    cout << "[INFO] smooth internal frames via L1" << endl;
    VectorXd fmat;
    convert_zyz_to_mat(abc, fmat);
    smoother->smoothL1(fmat);

    string ff_file = out_folder+string("/frames.txt");
    write_tet_ff(ff_file.c_str(), fmat.data(), fmat.size()/9);

  } else {
    cerr << "[Error] unsupported smooth term.\n";
    return __LINE__;
  }
  
  cout << "[Info] done\n";
  return 0;
}
