#include "lbfgs_solve.h"

#include <iostream>
#include <lbfgs.h>
#include <optimization.h>

#include "def.h"

using namespace std;
using namespace alglib;

namespace riemann {

static shared_ptr<Functional<double>> g_func;
static size_t g_count;

static lbfgsfloatval_t evaluate(void *instance, const lbfgsfloatval_t *x,
                                lbfgsfloatval_t *g, const int n,
                                const lbfgsfloatval_t step) {
  double value = 0;
  g_func->Val(x, &value);
  std::fill(g, g+n, 0);
  g_func->Gra(x, g);
  return value;
}

static int progress(void *instance,
                    const lbfgsfloatval_t *x,
                    const lbfgsfloatval_t *g,
                    const lbfgsfloatval_t fx,
                    const lbfgsfloatval_t xnorm,
                    const lbfgsfloatval_t gnorm,
                    const lbfgsfloatval_t step,
                    int n, int k, int ls) {
  if ( k % 100 == 0 )
    cout << "\t# iteration " << k << ", energy: " << fx << endl;
  //  printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
  //  printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
  //  printf("\n");
  return 0;
}

int lbfgs_solve(const shared_ptr<Functional<double>> &f,
                double *x, const size_t dim,
                const double epsf, const double epsx,
                const size_t maxiter) {
  if ( !f.get() ) {
    cerr << "[Error] null pointer to functional.\n";
    return __LINE__;
  }
  if ( dim != f->Nx() ) {
    cerr << "[Error] variable size not match.\n";
    return __LINE__;
  }
  
  g_func = f;
  g_count = 0;
  lbfgsfloatval_t fx;
  lbfgs_parameter_t param;
  lbfgs_parameter_init(&param);
  param.delta = epsf;
  param.max_iterations = maxiter;

  int ret = lbfgs(dim, x, &fx, evaluate, progress, NULL, &param);

  cout << "\t@LBFGS return value: " << ret << endl;
  return 0;
}

// static void function1_grad(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr) 
// {
//   const size_t dim = g_func->Nx();
//   const double *ptrx = x.getcontent();
//   double *ptrg = grad.getcontent();
  
//   func = 0;
//   g_func->Val(ptrx, &func);
//   //  cout << "\t@energy: " << func << endl;

//   std::fill(ptrg, ptrg+dim, 0);
//   g_func->Gra(ptrx, ptrg);
// }

// int lbfgs_solve(const shared_ptr<Functional<double>> &f,
//                 double *X, const size_t dim,
//                 const double EpsF, const double EpsX,
//                 const size_t maxiter) {
//   if ( !f.get() ) {
//     cerr << "[Error] null pointer to functional\n";
//     return __LINE__;
//   }
//   if ( f->Nx() != dim ) {
//     cerr << "[Error] variable size does not match\n";
//     return __LINE__;
//   }

//   g_func = f;
//   real_1d_array x;
//   x.setcontent(dim, X);
//   double epsg = 0.0000000001;
//   double epsf = EpsF;
//   double epsx = EpsX;
//   ae_int_t maxits = maxiter;
//   minlbfgsstate state;
//   minlbfgsreport rep;

//   minlbfgscreate(1, x, state);
//   minlbfgssetcond(state, epsg, epsf, epsx, maxits);
//   alglib::minlbfgsoptimize(state, function1_grad);
//   minlbfgsresults(state, x, rep);

//   const double *ptrx = x.getcontent();
//   std::copy(ptrx, ptrx+dim, X);

  
//   printf("\t# LBFGS RETURN: %d, ITERATIONS: %d\n", int(rep.terminationtype), int(rep.iterationscount));
//   return 0;
// }

}
