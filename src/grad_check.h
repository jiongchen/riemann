#ifndef GRAD_CHECK_H
#define GRAD_CHECK_H

// interface for user's customization
typedef void (*val_gra_func_t)(const int m, const int n, const double *x,
                               double *value, double *jacobian);

int numeric_grad_check(val_gra_func_t func, const int m, const int n, const double *x);

#endif
