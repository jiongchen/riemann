#include "cotmatrix.h"

#include "util.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

namespace surfparam {

void cotmatrix(const matrix<size_t> &cell,
               const matrix<double> &nods,
               const size_t dim,
               SparseMatrix<double> *L) {
    const size_t edge[3][2] = {{1, 2}, {2, 0}, {0, 1}};
    std::vector<Triplet<double>> trips;

    for (size_t i = 0; i < cell.size(2); ++i) {
        matrix<double> vert = nods(colon(), cell(colon(), i));
        matrix<double> half_cot_val(3);
        half_cot_val[0] = 0.5 * cal_cot_val(&vert(0, 1), &vert(0, 0), &vert(0, 2));
        half_cot_val[1] = 0.5 * cal_cot_val(&vert(0, 0), &vert(0, 1), &vert(0, 2));
        half_cot_val[2] = 0.5 * cal_cot_val(&vert(0, 0), &vert(0, 2), &vert(0, 1));

        for (size_t k = 0; k < 3; ++k) {
            size_t src = cell(edge[k][0], i);
            size_t des = cell(edge[k][1], i);
            runtime_dim_add_diag_block(dim, src, src, -half_cot_val[k], &trips);
            runtime_dim_add_diag_block(dim, des, des, -half_cot_val[k], &trips);
            runtime_dim_add_diag_block(dim, src, des, half_cot_val[k], &trips);
            runtime_dim_add_diag_block(dim, des, src, half_cot_val[k], &trips);
        }
    }
    const size_t lap_size = dim * nods.size(2);
    L->resize(lap_size, lap_size);
    L->reserve(trips.size());
    L->setFromTriplets(trips.begin(), trips.end());
}

}
