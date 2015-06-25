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
#include "src/nanoflann.hpp"

using namespace std;
using namespace Eigen;
using namespace zjucad::matrix;
using namespace surfparam;
using namespace geom_deform;
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
    geom_deform::quad_scalar_field_(&val, x, a, c);
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
    } catch (const boost::property_tree::ptree_error &e) {
        cerr << "Usage: " << endl;
        zjucad::show_usage_info(std::cerr, pt);
    } catch (const std::exception &e) {
        cerr << "# " << e.what() << endl;
    }
    return 0;
}
