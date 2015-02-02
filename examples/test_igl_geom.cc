#include <iostream>
#include <boost/filesystem.hpp>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/writeOBJ.h>
#include <igl/boundary_loop.h>

using namespace std;
using namespace Eigen;

MatrixXi tris;
MatrixXd nods;
MatrixXd uv;

int main(int argc, char *argv[])
{
    if ( argc != 2 ) {
        cerr << "# Usage: ./prog model.off\n";
        return __LINE__;
    }
    boost::filesystem::create_directory("./geom");
    igl::readOFF(argv[1], nods, tris);
    cout << "# tris size: " << tris.rows() << " " << tris.cols() << "\n";
    cout << "# nods size: " << nods.rows() << " " << nods.cols() << "\n";
    igl::writeOBJ("./geom/camelhead.obj", nods, tris);
    cout << "# done\n";
    return 0;
}
