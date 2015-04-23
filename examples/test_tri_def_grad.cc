#include <iostream>
#include <zjucad/matrix/io.h>
#include <jtflib/mesh/io.h>
#include <jtflib/mesh/util.h>
#include <boost/filesystem.hpp>
#include <hjlib/math/blas_lapack.h>
#include <zjucad/matrix/lapack.h>

#include "src/vtk.h"

using namespace std;
using namespace zjucad::matrix;

int main(int argc, char *argv[])
{
    boost::filesystem::create_directory("./tri_def_grad");

    matrix<size_t> tris = colon(0, 2);
    matrix<double> nods = rand(3, 3);
    jtf::mesh::save_obj("./tri_def_grad/origin.obj", tris, nods);

    matrix<double> normal;
    jtf::mesh::cal_face_normal(tris, nods, normal, true);
    matrix<size_t> line = colon(0, 1);
    matrix<double> line_nods(3, 2);
    line_nods(colon(), 0) = 1.0/3*nods*ones<double>(3, 1);
    line_nods(colon(), 1) = line_nods(colon(), 0)+normal;
    ofstream os("./tri_def_grad/normal.vtk");
    line2vtk(os, line_nods.begin(), line_nods.size(2), line.begin(), line.size(2));

    matrix<double> G(3, 3);
    G(colon(), colon(0, 1)) = nods(colon(), colon(1, 2))-nods(colon(), 0)*ones<double>(1, 2);
    G(colon(), 2) = normal;
    if ( inv(G) )
        cerr << "sigular matrix\n";

    matrix<double> def_nods = rand(3, 3);
    jtf::mesh::save_obj("./tri_def_grad/deform.obj", tris, def_nods);

    matrix<double> def_grad = zeros<double>(3, 3);
    def_grad(colon(), colon(0, 1)) = def_nods(colon(), colon(1, 2))-def_nods(colon(), 0)*ones<double>(1, 2);
    def_grad = temp(def_grad * G);

    matrix<double> dir0 = def_grad * (nods(colon(), 1)-nods(colon(), 0));
    matrix<double> dir1 = def_grad * (nods(colon(), 2)-nods(colon(), 0));
    matrix<size_t> edge_ele(2, 2);
    edge_ele(0, 0) = 0; edge_ele(0, 1) = 0;
    edge_ele(1, 0) = 1; edge_ele(1, 1) = 2;
    matrix<double> edge = zeros<double>(3, 3);
    edge(colon(), 1) = dir0;
    edge(colon(), 2) = dir1;
    ofstream ose("./tri_def_grad/edge.vtk");
    line2vtk(ose, edge.begin(), edge.size(2), edge_ele.begin(), edge_ele.size(2));

    cout << "done\n";
    return 0;
}
