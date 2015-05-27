#include <iostream>
#include <boost/filesystem.hpp>

#include "src/vec_field_deform.h"

using namespace std;
using namespace Eigen;

int main(int argc, char *argv[])
{
    if ( argc != 2 ) {
        cerr << "# Usage: " << argv[0] << " model.obj\n";
        return __LINE__;
    }
    boost::filesystem::create_directory("./vel_field_deform");

    geom_deform::vel_field_deform def;
    def.load_model(argv[1]);

//    // for sphere
//    Vector3d src(0, 1.0, 0);
//    Vector3d des(0, 1.2, 0);
//    const double ri = 0.01;
//    const double ro = 0.8;
//    def.translate_deform(src, des, ri, ro);

    // for bar
    Vector3d center(-6, 0, 0), n(1, 0, 0);
    const double ri = 12.05;
    const double ro = 14;
    def.twist_deform(center, ri, ro, n, 400);

    def.save_model("./vel_field_deform/deform.obj");
    cout << "done\n";
    return 0;
}
