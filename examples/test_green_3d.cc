#include <iostream>
#include <boost/filesystem.hpp>

#include "src/green_coord_deform.h"

using namespace std;
using namespace geom_deform;

int main(int argc, char *argv[])
{
    if ( argc != 3 ) {
        cerr << "usage: " << argv[0] << " model.obj cage.obj\n";
        return __LINE__;
    }
    boost::filesystem::create_directory("./green3d");

    green_deform_3d def;
    def.load_sample_points(argv[1]);
    def.load_cage(argv[2]);

    def.dump("./green/mesh.vtk");
    def.dump_cage("./green/cage.vtk");
    def.dump_normal("./green/cage_normal.vtk");

    cout << "done\n";
    return 0;
}
