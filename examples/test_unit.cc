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
#include "src/vec_field_deform.h"

using namespace std;
using namespace Eigen;
using namespace zjucad::matrix;
using namespace surfparam;
using namespace geom_deform;
using boost::property_tree::ptree;

int test_param_area(ptree &pt) {
    matrix<size_t> tris;
    matrix<double> nods, uv;
    jtf::mesh::load_obj("../../dat/lilium_param.obj", tris, nods);

    uv.resize(2, nods.size(2));
#pragma omp parallel for
    for (size_t i = 0; i < uv.size(2); ++i) {
        uv(0, i) = nods(0, i);
        uv(1, i) = nods(2, i);
    }

    shared_ptr<param_area> pa(new param_area(tris, nods));
    double area = 0;
    pa->Val(&uv[0], &area);
    cout << area << endl;

    cout << "# done\n";
    return 0;
}

//int test_quad_scalar_field(ptree &pt) {
//    double val = 0;
//    double a[3] = {3, 3, 0};
//    double x[3] = {1, 1, 0};
//    double c[3] = {2, 2, 0};
//    geom_deform::quad_scalar_field_(&val, x, a, c);
//    cout << "value: ";
//    cout << val << endl;
//    return 0;
//}

int test_vector_field(ptree &pt) {
    Vector3d c(0, 0, 0);
    Vector3d dir(1, 1, 1);
    const double ri = 1, ro = 2;
    vector_field vf(c, ri, ro, dir);

    Vector3d x(0, 1.5, 0);
    Vector3d vel = vf(x);
    cout << vel << endl;
    return 0;
}

int main(int argc, char *argv[])
{
    ptree pt;
    boost::filesystem::create_directory("./unitest");
    try {
        zjucad::read_cmdline(argc, argv, pt);
        CALL_SUB_PROG(test_param_area);
//        CALL_SUB_PROG(test_quad_scalar_field);
        CALL_SUB_PROG(test_vector_field);
    } catch (const boost::property_tree::ptree_error &e) {
        cerr << "Usage: " << endl;
        zjucad::show_usage_info(std::cerr, pt);
    } catch (const std::exception &e) {
        cerr << "# " << e.what() << endl;
    }
    return 0;
}
