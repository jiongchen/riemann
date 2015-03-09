#include <iostream>
#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

void cotmatrix(const matrix<size_t> &cell,
               const matrix<double> &nods,
               SparseMatrix<double> *L) {
    std::vector<Triplet<double>> trips;
    const size_t edge[3][2] = {{1, 2}, {2, 0}, {0, 1}};

    for (size_t i = 0; i < cell.size(2); ++i) {
        matrix<double> vert = nods(colon(), cell(colon(), i));

        const double area = [](const matrix<double> &V)->double {
            matrix<double> e = V(colon(), colon(1, 2)) - V(colon(), colon(0, 1));
            return 0.5 * norm(cross(e(colon(), 0), e(colon(), 1)));
        }(vert);

        matrix<double> len(3);
        for (size_t k = 0; k < 3; ++k) {
            len[i] = [&](const size_t id)->double {
                return norm(vert(colon(), cell(edge[id][0], i)) - vert(colon(), cell(edge[id][1], i)));
            }(k);
        }

        matrix<double> cotval(3);
//        cotval[0] = 1;
//        cotval[1] = 2;
//        cotval[2] = 3;

//        for (size_t k = 0; k < 3; ++k) {
//            size_t src = cell(edge[k][0], i);
//            size_t des = cell(edge[k][1], i);
//            add_diag_block();
//            add_diag_block();
//            add_diag_block();
//            add_diag_block();
//        }
    }

//    L->resize();
//    L->reserve();
//    L->setFromTriplets();
}

int main(int argc, char *argv[])
{
    for (size_t i = 0; i < 3; ++i) {
        [](const size_t t)->void {
            cout << t << endl;
        }(i);
    }
    return 0;
}
