// lake at rest with variable resolution
#include "saint-venant.h"

int main()
{
  origin (-0.5, -0.5);
  init_grid (16);
  run();
}

event init (i = 0)
{
  refine (level == 4 && x < -0.1 && y < -0.1);

  foreach() {
    zb[] = 0.2*exp(-100*(x*x + y*y));
    h[] = 1. - zb[];
  }
}

event logfile (i = 1)
{
  norm n = normf (u.x);
  fprintf (stderr, "# %.10f %.10f %.10f\n", n.avg, n.rms, n.max);
  foreach ()
    printf ("%g %g %g %g %g\n", x, y, h[], zb[], u.x[]);
}
