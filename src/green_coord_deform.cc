#include "green_coord_deform.h"

#include <zjucad/matrix/matrix.h>
#include <jtflib/mesh/io.h>

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

namespace geom_deform {

green_deform_2d::green_deform_2d() {}

int green_deform_2d::load_mesh(const char *file) {
    matrix<size_t> cell;
    matrix<double> nods;
    if ( jtf::mesh::load_obj(file, cell, nods) ) {
        cerr << "# info: input error\n";
        return __LINE__;
    }
    cell_.resize(cell.size(1), cell.size(2));
    std::copy(cell.begin(), cell.end(), cell_.data());

    ///> assign X
    X_.resize(2*nods.size(2));
#pragma omp parallel for
    for (size_t i = 0; i < nods_.size(2); ++i) {
        X_[2*i+0] = nods_(0, i);
        X_[2*i+1] = nods_(2, i);
    }
    return 0;
}

int green_deform_2d::load_cage(const char *file) {
    // parse manually
    ifstream is(file);
    size_t ele_num;
    is >> ele_num;
    cage_cell_.resize(2, ele_num);
    for (size_t i = 0; i < ele_num; ++i)
        is >> cage_cell_(0, i) >> cage_cell_(1, i);
    size_t nods_num;
    cage_nods_.resize(3, nods_num);
    for (size_t i = 0; i < nods_num; ++i);

    ///> assign N
    N_.resize(2*ele_num);
    calc_outward_normal();
    return 0;
}

int green_deform_2d::calc_outward_normal() {
#pragma omp parallel for
    for (size_t i = 0; i < cage_cell_.cols(); ++i) {
        Vector2d dir = Xcage_.segment<2>(2*cage_cell_(1, i))
                - Xcage_.segment<2>(2*cage_cell_(0, i));
        dir /= dir.norm();
    }
    return 0;
}

int green_deform_2d::calc_green_coords() {
    return 0;
}

int green_deform_2d::move_cage(const size_t id, const double *dx) {
    return 0;
}

int green_deform_2d::deform() {
    return 0;
}

int green_deform_2d::dump(const char *file) {
    return 0;
}

}
