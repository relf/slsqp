#ifndef SLSQP_H
#define SLSQP_H

#include <math.h>
#include <float.h>
// #include "nlopt.h"
// #include "nlopt-util.h"

// #ifdef __cplusplus
// extern "C"
// {
// #endif /* __cplusplus */

// nlopt_result nlopt_slsqp(unsigned n, nlopt_func f, void *f_data,
// 			 unsigned m, nlopt_constraint *fc,
// 			 unsigned p, nlopt_constraint *h,
// 			 const double *lb, const double *ub,
// 			 double *x, double *minf,
// 			 nlopt_stopping *stop);
// #ifdef __cplusplus
// }  /* extern "C" */
// #endif /* __cplusplus */


int nlopt_isnan(double x)
{
    return isnan(x);
}

int nlopt_isinf(double x)
{
    return (fabs(x) >= HUGE_VAL * 0.99)
        || isinf(x)
        ;
}

int nlopt_isfinite(double x)
{
    return (fabs(x) <= DBL_MAX)
        || isfinite(x)
        ;
}

typedef struct {
    double t, f0, h1, h2, h3, h4;
    int n1, n2, n3;
    double t0, gs;
    double tol;
    int line;
    double alpha;
    int iexact;
    int incons, ireset, itermx;
    double *x0;
} slsqpb_state;

void slsqp(int *m, int *meq, int *la, int *n,
		  double *x, const double *xl, const double *xu, double *f, 
		  double *c__, double *g, double *a, double *acc, 
		  int *iter, int *mode, double *w, int *l_w__, int *
		  jw, int *l_jw__, slsqpb_state *state);

#endif
