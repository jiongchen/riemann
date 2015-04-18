#include <iostream>
#include <boost/filesystem.hpp>

#include "src/green_coord_deform.h"

using namespace std;
using namespace geom_deform;

int main(int argc, char *argv[])
{
    if ( argc != 3 ) {
        cerr << "usage: " << argv[0] << " model.obj cage.2d\n";
        return __LINE__;
    }
    boost::filesystem::create_directory("./green");

    green_deform_2d def;
    def.load_sample_points(argv[1]);
    def.load_cage(argv[2]);
    def.calc_green_coords();
    {
        const double pos[2] = {-0.1, 1.1};
        def.move_cage(3, pos, false);
    }
    {
        const double pos[2] = {0.9, 1.1};
        def.move_cage(4, pos, false);
    }
    {
        const double pos[2] = {0.9, -0.1};
        def.move_cage(5, pos, false);
    }
    def.dump("./green/origin.vtk");
    def.deform();

    def.dump("./green/mesh.vtk");
    def.dump_cage("./green/cage.vtk");
    def.dump_normal("./green/cage_normal.vtk");

    cout << "done\n";
    return 0;
}
