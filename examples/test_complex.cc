#include <iostream>
#include <complex>
#include <cmath>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

using complexd=complex<double>;

int main(int argc, char *argv[])
{
  complexd ie(0, 1);
  cout << std::exp(ie*M_PI) << endl;
  Matrix<complexd, 4, 4> A;
  A = Matrix<complexd, 4, 4>::Random();
  cout << A << endl;
  cout << ie << endl;
  return 0;
}
