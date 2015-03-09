#include <igl/readOFF.h>
#include <igl/gaussian_curvature.h>
#include <igl/jet.h>

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  MatrixXd V;
  MatrixXi F;
  igl::readOFF(argv[1],V,F);

  VectorXd K;
  igl::gaussian_curvature(V,F,K);

//  cout << "guassian curvature: \n" << K << endl;
  int a = 0;
  [&a]()->void { a++;}();
  cout << a << endl;
  [](int &x)->void {x++;}(a);
  cout << a << endl;



  // Compute pseudocolor
  MatrixXd C;
  igl::jet(K,true,C);

  std::cout << "done\n";
  return 0;
}
