#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <zjucad/ptree/ptree.h>
#include <jtflib/mesh/mesh.h>
#include <jtflib/mesh/io.h>
#include <Eigen/Sparse>
#include <boost/filesystem.hpp>
#include <zjucad/matrix/io.h>
#include <igl/readOFF.h>
#include <igl/readDMAT.h>
#include <igl/grad.h>

#include "src/config.h"
#include "src/energy.h"
#include "src/vec_field_deform.h"
#include "src/nanoflann.hpp"
#include "src/cotmatrix.h"
#include "src/vtk.h"
#include "src/grad_operator.h"
#include "src/vert_local_frame.h"
#include "src/write_vtk.h"

using namespace std;
using namespace Eigen;
using namespace zjucad::matrix;
using namespace riemann;
using namespace riemann;
using namespace nanoflann;
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

int test_quad_scalar_field(ptree &pt) {
  double val = 0;
  double a[3] = {3, 3, 0};
  double x[3] = {1, 1, 0};
  double c[3] = {2, 2, 0};
  riemann::quad_scalar_field_(&val, x, a, c);
  cout << "value: ";
  cout << val << endl;
  return 0;
}

int test_vector_field(ptree &pt) {
  Vector3d c(0, 0, 0);
  Vector3d dir(1, 1, 1);
  const double ri = 1, ro = 2;
  vector_field vf(c, ri, ro, dir, "translate");

  Vector3d x(0, 1.5, 0);
  Vector3d vel = vf(x);
  cout << vel << endl;
  return 0;
}

const int SAMPLES_DIM = 15;

template <typename Der>
void generateRandomPointCloud(Eigen::MatrixBase<Der> &mat, const size_t N,const size_t dim, const typename Der::Scalar max_range = 10)
{
  std::cout << "Generating "<< N << " random points...";
  mat.resize(N,dim);
  for (size_t i=0;i<N;i++)
    for (size_t d=0;d<dim;d++)
      mat(i,d)= max_range * (rand() % 1000) / typename Der::Scalar(1000);
  std::cout << "done\n";
}

template <typename num_t>
void kdtree_demo(const size_t nSamples,const size_t dim)
{
  Eigen::Matrix<num_t,Dynamic,Dynamic>  mat(nSamples,dim);

  const num_t max_range = 20;

  // Generate points:
  generateRandomPointCloud(mat, nSamples,dim, max_range);

  //	cout << mat << endl;

  // Query point:
  std::vector<num_t> query_pt(dim);
  for (size_t d=0;d<dim;d++)
    query_pt[d] = max_range * (rand() % 1000) / num_t(1000);


  // ------------------------------------------------------------
  // construct a kd-tree index:
  //    Some of the different possibilities (uncomment just one)
  // ------------------------------------------------------------
  // Dimensionality set at run-time (default: L2)
  typedef KDTreeEigenMatrixAdaptor< Eigen::Matrix<num_t,Dynamic,Dynamic> >  my_kd_tree_t;

  // Dimensionality set at compile-time
  //	typedef KDTreeEigenMatrixAdaptor< Eigen::Matrix<num_t,Dynamic,Dynamic>, SAMPLES_DIM>  my_kd_tree_t;

  // Dimensionality set at compile-time: Explicit selection of the distance metric: L2
  //	typedef KDTreeEigenMatrixAdaptor< Eigen::Matrix<num_t,Dynamic,Dynamic>, SAMPLES_DIM,nanoflann::metric_L2>  my_kd_tree_t;

  // Dimensionality set at compile-time: Explicit selection of the distance metric: L2_simple
  //	typedef KDTreeEigenMatrixAdaptor< Eigen::Matrix<num_t,Dynamic,Dynamic>, SAMPLES_DIM,nanoflann::metric_L2_Simple>  my_kd_tree_t;

  // Dimensionality set at compile-time: Explicit selection of the distance metric: L1
  //	typedef KDTreeEigenMatrixAdaptor< Eigen::Matrix<num_t,Dynamic,Dynamic>, SAMPLES_DIM,nanoflann::metric_L1>  my_kd_tree_t;

  my_kd_tree_t   mat_index(dim /*dim*/, mat, 10 /* max leaf */ );
  mat_index.index->buildIndex();

  // do a knn search
  const size_t num_results = 3;
  vector<size_t>   ret_indexes(num_results);
  vector<num_t> out_dists_sqr(num_results);

  nanoflann::KNNResultSet<num_t> resultSet(num_results);

  resultSet.init(&ret_indexes[0], &out_dists_sqr[0] );
  mat_index.index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));

  std::cout << "knnSearch(nn="<<num_results<<"): \n";
  for (size_t i=0;i<num_results;i++)
    std::cout << "ret_index["<<i<<"]=" << ret_indexes[i] << " out_dist_sqr=" << out_dists_sqr[i] << endl;

  // do a radius search
  const num_t search_radius = static_cast<num_t>(250);
  std::vector<std::pair<long,num_t>> ret_matches;

  nanoflann::SearchParams params;
  //params.sorted = false;
  const size_t nMatches = mat_index.index->radiusSearch(&query_pt[0],search_radius, ret_matches, params);

  cout << "radiusSearch(): radius=" << search_radius << " -> " << nMatches << " matches\n";
  for (size_t i=0;i<nMatches;i++)
    cout << "idx["<< i << "]=" << ret_matches[i].first << " dist["<< i << "]=" << ret_matches[i].second << endl;
  cout << "\n";
}

int test_nanoflann(ptree &pt) {
  kdtree_demo<float>(1e3 /* samples */, SAMPLES_DIM /* dim */);
  return 0;
}

int test_matrix_extract(ptree &pt) {
  srand(time(NULL));

  matrix<double> A = rand(3, 4);
  cout << A << endl;

  matrix<size_t> p(2);
  p[0] = 1; p[1] = 3;

  matrix<double> B = rand(3, 2);
  cout << B << endl;

  A(colon(), p) = B;
  cout << A << endl;

  return 0;
}

int test_uniform_laplacian(ptree &pt) {
  matrix<size_t> cell;
  matrix<double> nods;
  jtf::mesh::load_obj("../../dat/plane_25.obj", cell, nods);
  nods += 0.1*rand(nods.size(1), nods.size(2));
  jtf::mesh::save_obj("./unitest/noise_model.obj", cell, nods);
  SparseMatrix<double> L;
  unimatrix(cell, nods, 3, &L, true);
  cout << L << endl;
  Map<VectorXd> X(nods.begin(), nods.size());
  int times = 1000;
  while ( times-- )
    X += 0.01*(L*X).eval();
  jtf::mesh::save_obj("./unitest/smoothed_model.obj", cell, nods);
  cout << "done\n";
  return 0;
}

int test_height_vector(ptree &pt) {
  srand(time(NULL));
  matrix<size_t> cell = colon(0, 2);
  matrix<double> nods = rand(3, 3);
  {
    ofstream os("./unitest/triangle.vtk");
    tri2vtk(os, nods.begin(), nods.size(2), cell.begin(), cell.size(2));
  }
  matrix<double> H(3, 3);
  riemann::calc_tri_height_vector<3>(&nods[0], &H[0]);
  matrix<size_t> line(2, 3);
  matrix<double> vert(3, 6);
  vert(colon(), colon(0, 2)) = nods;
  vert(colon(), colon(3, 5)) = nods-H;
  line(0, colon()) = colon(0, 2);
  line(1, colon()) = colon(3, 5);
  cout << line << endl;
  cout << vert << endl;
  {
    ofstream os("./unitest/height.vtk");
    line2vtk(os, &vert[0], vert.size(2), &line[0], line.size(2));
  }
  for (size_t j = 0; j < 3; ++j) {
    double len = dot(H(colon(), j), H(colon(), j));
    H(colon(), j) /= len;
  }
  cout << H*ones<double>(3, 1) << endl;
  cout << "done\n";
  return 0;
}

int test_grad_operator(ptree &pt) {
  matrix<size_t> tris;
  matrix<double> nods;
  jtf::mesh::load_obj("../../dat/half_sphere.obj", tris, nods);
  SparseMatrix<double> G;
  riemann::calc_grad_operator(tris, nods, &G);
  matd_t x = nods(0, colon());
  matd_t y = nods(1, colon());
  matd_t z = nods(2, colon());
  VectorXd grad_x = G*Map<VectorXd>(&y[0], y.size());
  itr_matrix<const double *> Gx(3, tris.size(2), grad_x.data());
  mati_t line(2, tris.size(2));
  matd_t vert(3, 2*tris.size(2));
  line(0, colon()) = colon(0, tris.size(2)-1);
  line(1, colon()) = colon(tris.size(2), 2*tris.size(2)-1);
#pragma omp parallel for
  for (size_t i = 0; i < tris.size(2); ++i) {
    vert(colon(), i) = nods(colon(), tris(colon(), i))*ones<double>(3, 1)/3.0;
    vert(colon(), i+tris.size(2)) = vert(colon(), i)+Gx(colon(), i);
  }
  jtf::mesh::save_obj("./unitest/test_model.obj", tris, nods);
  ofstream os("./unitest/grad.vtk");
  line2vtk(os, &vert[0], vert.size(2), &line[0], line.size(2));
  cout << "done\n";
  return 0;
}

int test_grad_operator2(ptree &pt) {
  Matrix<size_t, -1, -1> F;
  Matrix<double, -1, -1> V;
  igl::readOFF("../../dat/bunny.off", V, F);
  SparseMatrix<double> G0;
  igl::grad(V, F, G0);
  cout << G0.norm() << endl;
  cout << G0.blueNorm() << endl;
  cout << G0.nonZeros() << endl << endl;

  Matrix<size_t, -1, -1> FT = F.transpose();
  Matrix<double, -1, -1> VT = V.transpose();
  matrix<size_t> cell = itr_matrix<const size_t*>(FT.rows(), FT.cols(), FT.data());
  matrix<double> nods = itr_matrix<const double*>(VT.rows(), VT.cols(), VT.data());
  SparseMatrix<double> G1;
  riemann::calc_grad_operator(cell, nods, &G1);
  cout << G1.norm() << endl;
  cout << G1.blueNorm() << endl;
  cout << G1.nonZeros() << endl << endl;

  cout << "done\n";
  return 0;
}

int test_vert_local_frame(ptree &pt) {
  mati_t tris;
  matd_t nods;
  matd_t frame;
  jtf::mesh::load_obj("../../dat/beetle.obj", tris, nods);
  jtf::mesh::save_obj("./unitest/beetle.obj", tris, nods);
  riemann::calc_vert_local_frame(tris, nods, frame);
  frame *= 0.01;
  {
    matd_t x = frame(colon(0, 2), colon());
    riemann::draw_vert_direct_field("./unitest/tangent.vtk", nods.begin(), nods.size(2), x.begin());
  } {
    matd_t y = frame(colon(3, 5), colon());
    riemann::draw_vert_direct_field("./unitest/binormal.vtk", nods.begin(), nods.size(2), y.begin());
  } {
    matd_t z = frame(colon(6, 8), colon());
    riemann::draw_vert_direct_field("./unitest/normal.vtk", nods.begin(), nods.size(2), z.begin());
  }
  cout << "[info] done\n";
  return 0;
}

int main(int argc, char *argv[])
{
  ptree pt;
  boost::filesystem::create_directory("./unitest");
  try {
    zjucad::read_cmdline(argc, argv, pt);
    CALL_SUB_PROG(test_param_area);
    CALL_SUB_PROG(test_quad_scalar_field);
    CALL_SUB_PROG(test_vector_field);
    CALL_SUB_PROG(test_nanoflann);
    CALL_SUB_PROG(test_matrix_extract);
    CALL_SUB_PROG(test_uniform_laplacian);
    CALL_SUB_PROG(test_height_vector);
    CALL_SUB_PROG(test_grad_operator);
    CALL_SUB_PROG(test_grad_operator2);
    CALL_SUB_PROG(test_vert_local_frame);
  } catch (const boost::property_tree::ptree_error &e) {
    cerr << "Usage: " << endl;
    zjucad::show_usage_info(std::cerr, pt);
  } catch (const std::exception &e) {
    cerr << "# " << e.what() << endl;
  }
  return 0;
}
