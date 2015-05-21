#include <iostream>
#include <boost/filesystem.hpp>

#include "src/vec_field_deform.h"

#include <zjucad/linear_solver/linear_solver.h>

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
    Vector3d src(0, 1, 0);
    Vector3d des(0, 2, 0);
    const double ri = 0.05;
    const double ro = 0.3;
    def.translate(src, des, ri, ro);
    def.deform();
    def.save_model("./vel_field_deform/deform.obj");
    cout << "done\n";
    return 0;
}
