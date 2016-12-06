#include <iostream>
#include <memory>
#include <gtest/gtest.h>
#include <Eigen/Dense>
#include <IpIpoptApplication.hpp>
#include <IpSolveStatistics.hpp>
#include <limits>

#include "src/def.h"
#include "src/ipopt_solver.h"

using namespace std;
using namespace Eigen;
using namespace riemann;
using namespace Ipopt;

/**
 * argmin x0*x3*(x0+x1+x2)+x2
 * s.t.   x0*x1*x2*x3         >= 25
 *        x0^2+x1^2+x2^2+x3^2 =  40
 *        1 <= x0, x1, x2, x3 <= 5
 *
 * (1, 5, 5, 1) -> (...)
 */

class test_func : public Functional<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  test_func() {}
  size_t Nx() const {
    return 4;
  }
  int Val(const double *x, double *val) const {
    *val += x[0]*x[3]*(x[0]+x[1]+x[2])+x[2];
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    gra[0] += x[0]*x[3]+x[3]*(x[0]+x[1]+x[2]);
    gra[1] += x[0]*x[3];
    gra[2] += x[0]*x[3]+1;
    gra[3] += x[0]*(x[0]+x[1]+x[2]);
    return 0;
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    hes->push_back(TPL(0, 0, 2*x[3]));

    hes->push_back(TPL(0, 1, x[3]));
    hes->push_back(TPL(1, 0, x[3]));

    hes->push_back(TPL(0, 2, x[3]));
    hes->push_back(TPL(2, 0, x[3]));

    hes->push_back(TPL(0, 3, 2*x[0]+x[1]+x[2]));
    hes->push_back(TPL(3, 0, 2*x[0]+x[1]+x[2]));

    hes->push_back(TPL(1, 3, x[0]));
    hes->push_back(TPL(3, 1, x[0]));

    hes->push_back(TPL(2, 3, x[0]));
    hes->push_back(TPL(3, 2, x[0]));
    
    return 0;
  }
};

class test_cons: public Constraint<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  test_cons() {}
  size_t Nx() const {
    return 4;
  }
  size_t Nf() const {
    return 2;
  }
  int Val(const double *x, double *val) const {
    val[0] += x[0]*x[1]*x[2]*x[3];
    val[1] += x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3];
    return 0;
  }
  int Jac(const double *x, const size_t off, vector<Triplet<double>> *jac) const {    
    jac->push_back(TPL(off+0, 0, x[1]*x[2]*x[3]));
    jac->push_back(TPL(off+0, 1, x[0]*x[2]*x[3]));
    jac->push_back(TPL(off+0, 2, x[0]*x[1]*x[3]));
    jac->push_back(TPL(off+0, 3, x[0]*x[1]*x[2]));
    
    jac->push_back(TPL(off+1, 0, 2*x[0]));
    jac->push_back(TPL(off+1, 1, 2*x[1]));
    jac->push_back(TPL(off+1, 2, 2*x[2]));
    jac->push_back(TPL(off+1, 3, 2*x[3]));

    return 0;
  }
  int Hes(const double *x, const size_t off, vector<vector<Triplet<double>>> *hes) const {
    (*hes)[off+0].push_back(TPL(0, 1, x[2]*x[3]));
    (*hes)[off+0].push_back(TPL(0, 2, x[1]*x[3]));
    (*hes)[off+0].push_back(TPL(0, 3, x[1]*x[2]));
    (*hes)[off+0].push_back(TPL(1, 2, x[0]*x[3]));
    (*hes)[off+0].push_back(TPL(1, 3, x[0]*x[2]));
    (*hes)[off+0].push_back(TPL(2, 3, x[0]*x[1]));

    (*hes)[off+0].push_back(TPL(1, 0, x[2]*x[3]));
    (*hes)[off+0].push_back(TPL(2, 0, x[1]*x[3]));
    (*hes)[off+0].push_back(TPL(3, 0, x[1]*x[2]));
    (*hes)[off+0].push_back(TPL(2, 1, x[0]*x[3]));
    (*hes)[off+0].push_back(TPL(3, 1, x[0]*x[2]));
    (*hes)[off+0].push_back(TPL(3, 2, x[0]*x[1]));

    (*hes)[off+1].push_back(TPL(0, 0, 2));
    (*hes)[off+1].push_back(TPL(1, 1, 2));
    (*hes)[off+1].push_back(TPL(2, 2, 2));
    (*hes)[off+1].push_back(TPL(3, 3, 2));

    return 0;
  }
};

class test_nlp_problem : public ipopt_opt_framework
{
public:
  test_nlp_problem(const shared_ptr<Functional<double>> &obj,
                   const shared_ptr<Constraint<double>> &con,
                   double *x0)
      : ipopt_opt_framework(obj, con, x0) {}
  bool get_bounds_info(Index n, Number *x_l, Number *x_u,
                       Index m, Number *g_l, Number *g_u) {
    for (size_t i = 0; i < n; ++i) {
      x_l[i] = 1.0;
      x_u[i] = 5.0;
    }

    g_l[0] = 25;
    g_u[0] = 2e19;

    g_l[1] = 40.0;
    g_u[1] = 40.0;
    
    return true;
  }
};
  
TEST(ipopt_wrapper, with_inequality)
{
  shared_ptr<Functional<double>> obj = make_shared<test_func>();

  vector<shared_ptr<Constraint<double>>> cbuf(1);
  cbuf[0] = make_shared<test_cons>();
  shared_ptr<Constraint<double>> cons;
  try {
    cons = make_shared<constraint_t<double>>(cbuf);
  } catch (...) {
    cerr << "# exceptions" << endl;
    exit(EXIT_FAILURE);
  }
  
  // solve
  double x0[4] = {1, 5, 5, 1};
  
  SmartPtr<TNLP> mynlp = new test_nlp_problem(obj, cons, x0);
  
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
  app->RethrowNonIpoptException(true);

  app->Options()->SetNumericValue("tol", 1e-7);
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetStringValue("output_file", "ipopt.out");
      
  ApplicationReturnStatus status;
  status = app->Initialize();
  
  if (status != Solve_Succeeded) {
    std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
    exit(1);
  }

  status = app->OptimizeTNLP(mynlp);
  
  if (status == Solve_Succeeded) {
    std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
  } else {
    std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
  }

  cout << "# solution: " << endl;
  for (size_t i = 0; i < 4; ++i)
    cout << x0[i] << endl;
  
  EXPECT_EQ((int)status, 0);
}

class test_func2 : public Functional<double>
{
public:
  typedef Eigen::Triplet<double> TPL;
  test_func2() {}
  size_t Nx() const {
    return 4;
  }
  int Val(const double *x, double *val) const {
    *val = pow(x[0]-3, 2)+pow(x[1]+1, 2)+pow(x[2]-M_PI, 2)+pow(x[3]+exp(1), 2);
    return 0;
  }
  int Gra(const double *x, double *gra) const {
    gra[0] = 2*(x[0]-3);
    gra[1] = 2*(x[1]+1);
    gra[2] = 2*(x[2]-M_PI);
    gra[3] = 2*(x[3]+exp(1));
    return 0;
  }
  int Hes(const double *x, vector<Triplet<double>> *hes) const {
    hes->push_back(TPL(0, 0, 2));
    hes->push_back(TPL(1, 1, 2));
    hes->push_back(TPL(2, 2, 2));
    hes->push_back(TPL(3, 3, 2));
    return 0;
  }
};

TEST(ipopt_wrapper, unconstrained)
{
  shared_ptr<Functional<double>> obj = make_shared<test_func2>();
  shared_ptr<Constraint<double>> null_con;
  
  // solve
  double x0[4] = {1, 5, 5, 1};
  
  SmartPtr<TNLP> mynlp = new ipopt_opt_framework(obj, null_con, x0);
  
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
  app->RethrowNonIpoptException(true);

  app->Options()->SetNumericValue("tol", 1e-7);
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetStringValue("output_file", "ipopt.out");
      
  ApplicationReturnStatus status;
  status = app->Initialize();
  
  if (status != Solve_Succeeded) {
    std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
    exit(1);
  }

  status = app->OptimizeTNLP(mynlp);
  
  if (status == Solve_Succeeded) {
    std::cout << std::endl << std::endl << "*** The problem solved!" << std::endl;
  } else {
    std::cout << std::endl << std::endl << "*** The problem FAILED!" << std::endl;
  }

  cout << "# solution: " << endl;
  for (size_t i = 0; i < 4; ++i)
    cout << x0[i] << endl;
  
  EXPECT_EQ((int)status, 0);
}


TEST(ipopt_wrapper, int32andptrdiff) {
  EXPECT_EQ(sizeof(int64_t), sizeof(ptrdiff_t));
  EXPECT_EQ(std::numeric_limits<int64_t>::min(), std::numeric_limits<ptrdiff_t>::min());
  EXPECT_EQ(std::numeric_limits<int64_t>::max(), std::numeric_limits<ptrdiff_t>::max());
}
