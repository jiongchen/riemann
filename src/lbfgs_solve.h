#ifndef LBFGS_SOLVE_H
#define LBFGS_SOLVE_H

#include <memory>
#include <optimization.h>

namespace riemann {

template <typename T>
class Functional;

/*************************************************************************
This function sets stopping conditions for L-BFGS optimization algorithm.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    EpsG    -   >=0
                The  subroutine  finishes  its  work   if   the  condition
                |v|<EpsG is satisfied, where:
                * |.| means Euclidian norm
                * v - scaled gradient vector, v[i]=g[i]*s[i]
                * g - gradient
                * s - scaling coefficients set by MinLBFGSSetScale()
    EpsF    -   >=0
                The  subroutine  finishes  its work if on k+1-th iteration
                the  condition  |F(k+1)-F(k)|<=EpsF*max{|F(k)|,|F(k+1)|,1}
                is satisfied.
    EpsX    -   >=0
                The subroutine finishes its work if  on  k+1-th  iteration
                the condition |v|<=EpsX is fulfilled, where:
                * |.| means Euclidian norm
                * v - scaled step vector, v[i]=dx[i]/s[i]
                * dx - ste pvector, dx=X(k+1)-X(k)
                * s - scaling coefficients set by MinLBFGSSetScale()
    MaxIts  -   maximum number of iterations. If MaxIts=0, the  number  of
                iterations is unlimited.

Passing EpsG=0, EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to
automatic stopping criterion selection (small EpsX).

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/

/*************************************************************************
 2922 L-BFGS algorithm results
 2923 
 2924 INPUT PARAMETERS:
 2925     State   -   algorithm state
 2926 
 2927 OUTPUT PARAMETERS:
 2928     X       -   array[0..N-1], solution
 2929     Rep     -   optimization report:
 2930                 * Rep.TerminationType completetion code:
 2931                     * -7    gradient verification failed.
 2932                             See MinLBFGSSetGradientCheck() for more information.
 2933                     * -2    rounding errors prevent further improvement.
 2934                             X contains best point found.
 2935                     * -1    incorrect parameters were specified
 2936                     *  1    relative function improvement is no more than
 2937                             EpsF.
 2938                     *  2    relative step is no more than EpsX.
 2939                     *  4    gradient norm is no more than EpsG
 2940                     *  5    MaxIts steps was taken
 2941                     *  7    stopping conditions are too stringent,
 2942                             further improvement is impossible
 2943                 * Rep.IterationsCount contains iterations count
 2944                 * NFEV countains number of function calculations
 2945 
 2946   -- ALGLIB --
 2947      Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/

int lbfgs_solve(const std::shared_ptr<Functional<double>> &f,
                double *X, const size_t dim,
                const double EpsF = 0, const double EpsX = 0,
                const size_t maxiter = 0);

using alglib::real_1d_array;

typedef void (*alglib_callback_t)(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr);

int lbfgs_solve(alglib_callback_t cb,
                double *X, const size_t dim,
                const double EpsF = 0, const double EpsX = 0,
                const size_t maxiter = 0);
}

#endif
