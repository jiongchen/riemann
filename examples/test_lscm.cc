#include <iostream>
#include <boost/filesystem.hpp>

#include "igl/readOBJ.h"
#include "igl/readOFF.h"
#include "igl/writeOBJ.h"
#include "igl/boundary_loop.h"
#include "igl/lscm.h"

using namespace std;
using namespace Eigen;

MatrixXd V;
MatrixXi F;
MatrixXd uv;
MatrixXd UV;

int main(int argc, char *argv[])
{
    if ( argc != 2 ) {
        cerr << "# Uasge: ./prog model.off\n";
        return __LINE__;
    }
    igl::readOBJ(argv[1], V, F);

    VectorXi bnd, b(2, 1);
    igl::boundary_loop(F, bnd);
    b(0) = bnd(0);
    b(1) = bnd(round(bnd.size() / 2));
    cout << "[INFO] boundary vert: " << b(0) << " " << b(1) << endl;
    MatrixXd bc(2, 2);
    bc << 0, 0, 1, 0;

    // Least Square Conformal Mapping parametrization
    igl::lscm(V, F, b, bc, uv);
    uv *= 5;

    boost::filesystem::create_directory("./lscm");
    UV.resize(uv.rows(), 3);
    UV.setZero();
    UV.block(0, 0, uv.rows(), uv.cols()) = uv;
    igl::writeOBJ("./lscm/param.obj", UV, F);

    cout << "# done\n";
    return 0;
}
