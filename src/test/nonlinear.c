#include <gsl/gsl_integration.h>
#pragma autolink -lgsl -lgslcblas
#include "grid/cartesian.h"
#include "atmosphere.h"

int main() {
  // coordinates of lower-left corner
  origin (-0.5, -0.5);
  // number of grid points
  init_grid (128);
  // size of the box
  size (1.);
  // acceleration of gravity
  G = 1.;
  // Coriolis parameter
  F0 = 0.;
  // Viscosity
  //  NU = 7e-5;
  NU = 0.;
  // CFL number
  //  CFL = 0.5;
  CFL = 0.25;
  run();
}

/* ---------------- Initial conditions ------------------- */

#define H0 1.
#define FROUDE 0.1

double vtheta (double r) {
  return FROUDE*(r < 0.4)*(1. + cos((r - 0.2)/0.2*M_PI))/2.;
}

double h0p (double r, void * p) {
  double vt = vtheta(r);
  return vt*(F0 + vt/r)/G;
}

double h0 (double r) {
  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  double result, error;
  gsl_function F;
  F.function = &h0p;
  gsl_integration_qags (&F, 0, r, 0, 1e-7, 1000, w, &result, &error);
  gsl_integration_workspace_free (w);
  return result;
}

scalar h1[];

event init (i = 0)
{
  foreach()
    h1[] = h[] = (H0 + h0(sqrt (x*x + y*y)));
  foreach_face(x)
    u.x[] = - vtheta(sqrt (x*x + y*y))*y/sqrt (x*x + y*y);
  foreach_face(y)
    u.y[] =   vtheta(sqrt (x*x + y*y))*x/sqrt (x*x + y*y);
}

/* ------------------ Output helper functions --------------- */

scalar e[];

double error()
{
  double max = 0.;
  foreach(reduction(max:max)) {
    e[] = fabs (h1[]  - h[]);
    if (e[] > max) max = e[];
  }
  return max;
}

double energy()
{
  double se = 0.;
  foreach(reduction(+:se)) {
    double ke = (sq(u.x[] + u.x[1,0]) + sq(u.y[] + u.y[0,1]))/8.;
    se += (h[]*ke + G*h[]*(h[]/2. + b[]))*sq(Delta);
  }
  return se;
}

event logfile (i += 10; t <= 5.)
  fprintf (stderr, "t: %g %g %g\n", t, error(), energy());
