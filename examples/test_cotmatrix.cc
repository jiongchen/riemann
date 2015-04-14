#include <iostream>
#include <zjucad/matrix/matrix.h>
#include <Eigen/Sparse>
#include <unsupported/Eigen/KroneckerProduct>

#include "igl/cotmatrix.h"
#include "igl/readOBJ.h"

using namespace std;
using namespace zjucad::matrix;
using namespace Eigen;

int main(int argc, char *argv[])
{
    MatrixXi F;
    MatrixXd V;
    igl::readOBJ("../../dat/sphere.obj", V, F);

    SparseMatrix<double> L;
    igl::cotmatrix(V, F, L);

    SparseMatrix<double> Lap = kroneckerProduct(L, Matrix3d::Identity());
    cout << Lap.coeff(200, 200) << endl;
    cout << "done\n";
    return 0;
}
