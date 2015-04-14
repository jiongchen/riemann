#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <zjucad/ptree/ptree.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <Eigen/Sparse>
#include <boost/filesystem.hpp>
#include <zjucad/matrix/io.h>

#include "src/config.h"
#include "src/energy.h"

using namespace std;
using namespace Eigen;
using namespace zjucad::matrix;
using namespace surfparam;
using boost::property_tree::ptree;

int test_param_area(ptree &pt) {
    matrix<size_t> tris;
    matrix<double> nods, uv;
    jtf::mesh::load_obj("../../dat/plane_with_hole.obj", tris, nods);

    uv.resize(2, nods.size(2));
#pragma omp parallel for
    for (size_t i = 0; i < uv.size(2); ++i) {
        uv(0, i) = nods(0, i);
        uv(1, i) = nods(2, i);
    }

    shared_ptr<param_area> pa(new param_area(tris, nods, 1.0));
    double area = 0;
    pa->Val(&uv[0], &area);
    cout << area << endl;

    cout << "# done\n";
    return 0;
}

int main(int argc, char *argv[])
{
    ptree pt;
    boost::filesystem::create_directory("./unitest");
    try {
        zjucad::read_cmdline(argc, argv, pt);
        CALL_SUB_PROG(test_param_area);
    } catch (const boost::property_tree::ptree_error &e) {
        cerr << "Usage: " << endl;
        zjucad::show_usage_info(std::cerr, pt);
    } catch (const std::exception &e) {
        cerr << "# " << e.what() << endl;
    }
    return 0;
}
