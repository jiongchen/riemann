#include <iostream>
#include <fstream>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <boost/program_options.hpp>

#include "src/volume_frame.h"
#include "src/vtk.h"
#include "src/lbfgs_solve.h"

using namespace std;
using namespace riemann;
using namespace Eigen;
namespace po=boost::program_options;

class test_func : public Functional<double>
{
public:
  size_t Nx() const {
    return 2;
  }
  int Val(const double *x, double *val) const {
    *val += 100*pow(x[0]+3,4) + pow(x[1]-3,4);
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    gra[0] += 400*pow(x[0]+3,3);
    gra[1] += 4*pow(x[1]-3,3);
    return 0;
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    return __LINE__;
  }
};

static int simple_test_on_lbfgs()
{
  shared_ptr<Functional<double>> F = make_shared<test_func>();
  double x[2] = {0};
  
  const double epsf = 0, epsx = 0;
  const size_t maxits = 0;
  lbfgs_solve(F, x, 2, epsf, epsx, maxits);

  cout << x[0] << " " << x[1] << endl;
  return 0;
}

int main(int argc, char *argv[])
{
  po::options_description desc("Available options"); {
    desc.add_options()
        ("help,h", "produce help message")
        ("mesh,i", po::value<string>(), "input mesh")
        ("output_folder,o", po::value<string>(), "output folder")
        ;
  }
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if ( vm.count("help") ) {
    cout << desc << endl;
    return __LINE__;
  }

  mati_t tets; matd_t nods;
  jtf::mesh::tet_mesh_read_from_zjumat(vm["mesh"].as<string>().c_str(), &nods, &tets);

  shared_ptr<cross_frame_opt> frame_opt(cross_frame_opt::create(tets, nods));

  VectorXd Fs = VectorXd::Zero(9*nods.size(2));
  frame_opt->solve_smooth_sh_coeffs(Fs);

  VectorXd abc = VectorXd::Zero(3*nods.size(2));
  frame_opt->solve_initial_frames(Fs, abc);
  frame_opt->optimize_frames(abc);


  cout << "[Info] done\n";
  return 0;
}
