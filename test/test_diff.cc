#include <iostream>
#include <memory>
#include <gtest/gtest.h>
#include <Eigen/Dense>

#include "src/grad_check.h"

using namespace std;
using namespace Eigen;

class test_func
{
public:
  int dim() const {
    return 3;
  }
  void val(const double *x, double *f) const {
    *f = std::pow(x[0]-1, 3)+std::pow(x[1]-2, 2)+std::sin(x[2]);
  }
  void gra(const double *x, double *g) const {
    g[0] = 3*std::pow(x[0]-1, 2);
    g[1] = 2*(x[1]-2);
    g[2] = std::cos(x[2]);
  }
};

static shared_ptr<test_func> gf;

static void evaluate(const int m, const int n, const double *x,
                     double *value, double *jacobian) {
  gf->val(x, value);
  gf->gra(x, jacobian);
}

TEST(diff_test, gradient) {
  gf = make_shared<test_func>();
  int m = 1;
  int n = gf->dim();

  srand(time(NULL));
  VectorXd x = VectorXd::Random(n);

  int rtn = numeric_grad_check(evaluate, m, n, x.data());

  EXPECT_EQ(rtn, 0);
}
