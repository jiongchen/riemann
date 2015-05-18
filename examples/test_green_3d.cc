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
    def.calc_green_coords();

//    /// for sphere
//    const double disp[3] = {0.0, -2.0, 0.0};
//    def.move_cage(4, disp, true);
//    def.move_cage(5, disp, true);
//    def.move_cage(6, disp, true);
//    def.move_cage(7, disp, true);

    /// for bar
    const double pos4[3] = {-1.031485, 4.966156, 1.048004};
    const double pos5[3] = {-1.031485, 4.966156, -1.046975};
    const double pos6[3] = {1.063854, 4.966156, -1.046975};
    const double pos7[3] = {1.063854, 4.966156, 1.048004};
    def.move_cage(4, pos5, false);
    def.move_cage(5, pos6, false);
    def.move_cage(6, pos7, false);
    def.move_cage(7, pos4, false);

    def.deform();

    def.dump("./green3d/mesh.vtk");
    def.dump_cage("./green3d/cage.vtk");
    def.dump_normal("./green3d/cage_normal.vtk");

    cout << "done\n";
    return 0;
}
