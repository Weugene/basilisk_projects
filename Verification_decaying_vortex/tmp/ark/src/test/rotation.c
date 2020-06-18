#include "grid/cartesian.h"
#include "advection.h"

scalar f[];
scalar * tracers = {f};

u.n[left]   = 0.;
u.n[right]  = 0.;
u.n[top]    = 0.;
u.n[bottom] = 0.;

int main()
{
  // coordinates of lower-left corner
  origin (-0.5, -0.5);
  // maximum timestep
  DT = .1;
  // CFL number
  CFL = 0.8;
  for (N = 64; N <= 256; N *= 2)
    run ();
}

double bump (double x, double y)
{
  double r2 = x*x + y*y; 
  double coeff = 20. + 20000.*r2*r2*r2*r2;
  return (1. + cos(20.*x)*cos(20.*y))*exp(-coeff*r2)/2.;
}

event init (i = 0)
{
  foreach()
    f[] = bump(x,y);
}

#define end 0.785398

event velocity (i++) {
  trash ({u});
  foreach_face(x) u.x[] = -8.*y;
  foreach_face(y) u.y[] =  8.*x;
  boundary ((scalar *){u});
}

event logfile (t = {0,end}) {
  stats s = statsf (f);
  fprintf (stderr, "# %f %.12f %g %g\n", t, s.sum, s.min, s.max);
}

event field (t = end) {
  scalar e[];
  foreach()
    e[] = f[] - bump(x,y);
  norm n = normf (e);
  fprintf (stderr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);

  if (N == 256)
    output_field ({e}, stdout, N);
}
