#include "grid/cartesian.h"
#include "advection.h"

scalar f[];
scalar * tracers = {f};

int main()
{
  // coordinates of lower-left corner
  origin (-0.5, -0.5);
  // maximum timestep
  DT = .1;
  // CFL number
  CFL = 0.8;
  for (N = 64; N <= 256; N *= 2)
    run();
}

#define bump(x,y) (exp(-100.*(sq(x + 0.2) + sq(y + .236338))))

event init (i = 0)
{
  foreach()
    f[] = bump(x,y);
}

event velocity (i++) {
  vertex scalar psi[];
  foreach_vertex()
    psi[] = - 1.5*sin(2.*pi*t/5.)*sin((x + 0.5)*pi)*sin((y + 0.5)*pi)/pi;
  trash ({u});
  struct { double x, y; } f = {-1.,1.};
  foreach_face()
    u.x[] = f.x*(psi[0,1] - psi[])/Delta;
  boundary ((scalar *){u});
}

event logfile (t = {0,5}) {
  stats s = statsf (f);
  fprintf (stderr, "# %f %.12f %g %g\n", t, s.sum, s.min, s.max);
}

event field (t = 5) {
  scalar e[];
  foreach()
    e[] = f[] - bump(x,y);
  norm n = normf (e);
  fprintf (stderr, "%d %g %g %g\n", N, n.avg, n.rms, n.max);
  
  if (N == 256)
    output_field ({e}, stdout, N);
}
