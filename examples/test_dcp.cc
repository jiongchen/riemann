#include <iostream>
#include <boost/filesystem.hpp>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/mesh.h>

#include "src/lscm.h"

using namespace std;
using namespace zjucad::matrix;
using namespace surfparam;

int main(int argc, char *argv[])
{
    if ( argc != 2 ) {
        cerr << "usage: " << argv[0] << " model.obj\n";
        return __LINE__;
    }
    boost::filesystem::create_directory("./dcp");

    matrix<size_t> tris, bnd_e;
    matrix<double> nods;
    if ( jtf::mesh::load_obj(argv[1], tris, nods) ) {
        cerr << "# no input file\n";
        return __LINE__;
    }
    jtf::mesh::reorder_face(tris);
    jtf::mesh::save_obj("./dcp/orgin.obj", tris, nods);

    shared_ptr<jtf::mesh::edge2cell_adjacent> e2c(jtf::mesh::edge2cell_adjacent::create(tris, false));
    jtf::mesh::get_boundary_edge(*e2c, bnd_e);
    sort(bnd_e.begin(), bnd_e.end());

    lscm_param param(tris, nods);
    const double pos[4] = {0, 0, 1, 0};
    param.set_fixed_bnd_vert(bnd_e[0], &pos[0]);
    param.set_fixed_bnd_vert(bnd_e[bnd_e.size()/2], &pos[2]);
    cout << "# fixed point: " << bnd_e[0] << " " << bnd_e[bnd_e.size()/2] << endl;

    param.apply();

    matrix<size_t> param_tris;
    matrix<double> param_nods;
    param.get_param_mesh(&param_tris, &param_nods);
    param_nods *= 5.0;
    jtf::mesh::save_obj("./dcp/param.obj", param_tris, param_nods);

    cout << "# done\n";
    return 0;
}
