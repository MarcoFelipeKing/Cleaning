/* file decay_ODE.c */
#include <R.h>
static double parms[4];
#define C parms[0]
#define d parms[1]
#define g parms[2]
#define r parms[3]

/* initializer  */
void initmod(void (* odeparms)(int *, double *))
{
  int N=4;
  odeparms(&N, parms);
}

/* Derivatives and 1 output variable */
void derivs (int *neq, double *t, double *y, double *ydot,
             double *yout, int *ip)
{
  // if (ip[0] <1) error("nout should be at least 1");
  ydot[0] = y[0]*(1 - y[0] / C) * r - d * exp(-g * *t) * y[0];
}
/* END file decayODE.c */
