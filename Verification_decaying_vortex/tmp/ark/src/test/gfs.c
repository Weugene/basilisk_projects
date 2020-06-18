#include "utils.h"

int main()
{
  int depth = 6;
  origin (-0.5, -0.5, -0.5);
  init_grid (8);
  refine (level < depth - 2 || level <= depth*(1. - sqrt(x*x + y*y + z*z)));
  
  scalar a[];
  vector u[];
  foreach()
    u.x[] = u.y[] = a[] = sin(2.*pi*x)*cos(2.*pi*y);
  boundary ({a,u});

  FILE * fp = fopen ("gfsi.gfs", "w");
  output_gfs (fp);
  fclose (fp);
}
