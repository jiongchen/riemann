#include <iostream>
#include <boost/filesystem.hpp>
#include <jtflib/mesh/io.h>
#include <Eigen/Geometry>

#include "src/arap_deform.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

int main(int argc, char *argv[])
{
    if ( argc != 2 ) {
        cerr << "Uasge: " << argv[0] << " model.obj\n";
        return __LINE__;
    }
    boost::filesystem::create_directory("./arap");
    matrix<size_t> tris;
    matrix<double> nods;
    jtf::mesh::load_obj(argv[1], tris, nods);
    jtf::mesh::save_obj("./arap/origin.obj", tris, nods);

    core::arap_deform def(tris, nods);

#define TEST1 1
#if TEST1
    vector<size_t> idx{0, 1, 2, 3, 4, 5, 6, 7};
    Matrix3d rot;
    rot = AngleAxisd(M_PI/2, Vector3d::UnitX());
    matrix<double> R(3, 3);
    std::copy(rot.data(), rot.data()+9, R.begin());
    nods(colon(), 0) = temp(R*nods(colon(), 0));
    nods(colon(), 1) = temp(R*nods(colon(), 1));
    nods(colon(), 5) = temp(R*nods(colon(), 5));
    nods(colon(), 4) = temp(R*nods(colon(), 4));
# endif
#if TEST2
    vector<size_t> idx{15, 391};
    matrix<double> disp(3);
    disp[0] = disp[2] = 0.0;
    disp[1] = 1.0;
    nods(colon(), 391) += disp;
#endif
#define TEST3 0
#if TEST3
    vector<size_t> idx{0, 1, 4, 9, 14, 2, 3, 5, 15, 10};
    matrix<double> disp(3);
    disp[0] = 0.0;
    disp[1] = 1.0;
    disp[2] = -1.0;
    nods(colon(), 0) += disp;
    nods(colon(), 1) += disp;
    nods(colon(), 4) += disp;
    nods(colon(), 9) += disp;
    nods(colon(), 14) += disp;
#endif
#define TEST4 0
#if TEST4
    vector<size_t> idx{2, 31, 15, 47, 5, 42, 10, 26, 3,
                       0, 30, 14, 46, 4, 41, 9, 25, 1};
    matrix<double> disp(3);
    disp[0] = disp[1] = 0.0;
    disp[2] = -1.5;
    for (size_t i = 9; i < idx.size(); ++i)
        nods(colon(), idx[i]) += disp;
#endif

    def.pre_compute(idx);
    def.deformation(&nods[0]);

    jtf::mesh::save_obj("./arap/deform.obj", tris, nods);
    cout << "# done\n";
    return 0;
}
