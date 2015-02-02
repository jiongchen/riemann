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

int main(int argc, char *argv[])
{
    if ( argc != 2 ) {
        cerr << "# Usage: ./prog model.off\n";
        return __LINE__;
    }
    boost::filesystem::create_directory("./geom");
    igl::readOFF(argv[1], nods, tris);
    igl::writeOBJ("./geom/camelhead.obj", nods, tris);

    cout << "# done\n";
    return 0;
}
