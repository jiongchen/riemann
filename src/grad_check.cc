#include "grad_check.h"

#include <iostream>
#include <minpack.h>
#include <cstdlib>

using namespace std;

struct chkder_prb {
  int m, n;
  double *x;     // length n
  double *fvec;  // length m
  double *fjac;  // length m by n
  int ldfjac;    // >= m
  double *xp;    // length n
  double *fvecp; // length m
  int mode;
  double *err;   // length m
};

int numeric_grad_check(val_gra_func_t func, const int m, const int n, const double *x) {
  chkder_prb prb;
  prb.m = m;
  prb.n = n;
  prb.x = x;

  // allocate memory
  prb.fvec  = (double *)malloc(m*sizeof(double));
  prb.fjac  = (double *)malloc(m*n*sizeof(double));
  prb.ldfjac = m;
  prb.xp    = (double *)malloc(n*sizeof(double));
  prb.fvecp = (double *)malloc(m*sizeof(double));
  prb.err   = (double *)malloc(m*sizeof(double));
  
  prb.mode = 1;
  chkder_(&prb.m, &prb.n, prb.x, prb.fvec, prb.fjac, &prb.ldfjac,
          prb.xp, prb.fvecp, &prb.mode, prb.err);

  // evalute fvecp and fjac;
  func(prb.m, prb.n, prb.xp, prb.fvecp, prb.fjac);
  func(prb.m, prb.n, prb.x, prb.fvec, prb.fjac);
  
  prb.mode = 2;
  chkder_(&prb.m, &prb.n, prb.x, prb.fvec, prb.fjac, &prb.ldfjac,
          prb.xp, prb.fvecp, &prb.mode, prb.err);

  int flag = 0;
  for (int i = 0; i < m; ++i) {
    if ( prb.err[i] < 0.5 ) {
      cout << "## f[" << i << "]=" << prb.err[i] << " is probably incorrect." << endl;
      flag = 1;
    }
  }
  
  free(prb.fvec);
  free(prb.fjac);
  free(prb.xp);
  free(prb.fvecp);
  free(prb.err);

  return flag;
}
