#include <iostream>
#include <boost/filesystem.hpp>

#include "igl/readOBJ.h"
#include "igl/readOFF.h"
#include "igl/writeOBJ.h"
#include "igl/boundary_loop.h"
#include "igl/map_vertices_to_circle.h"
#include "igl/harmonic.h"

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
    igl::readOFF(argv[1], V, F);

    VectorXi bnd;
    igl::boundary_loop(V, F, bnd);

    MatrixXd bnd_uv;
    igl::map_vertices_to_circle(V, bnd, bnd_uv);

    igl::harmonic(V, F, bnd, bnd_uv, 1, uv);
    uv *= 5;

    boost::filesystem::create_directory("./harmonic");
    UV.resize(uv.rows(), 3);
    UV.setZero();
    UV.block(0, 0, uv.rows(), uv.cols()) = uv;
    igl::writeOBJ("./harmonic/param.obj", UV, F);

    cout << "# done\n";
    return 0;
}
