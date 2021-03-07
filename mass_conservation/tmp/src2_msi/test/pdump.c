// generates the dump file read by restore.c

#include "utils.h"

int main()
{
  int depth = 6;
  origin (-0.5, -0.5, -0.5);

#if TREE
  init_grid (1);
  refine (level < depth - 2 || level <= depth*(1. - sqrt(x*x + y*y + z*z)));
#else
  init_grid (1 << depth);
#endif
  
  scalar s[];
  foreach()
    s[] = sin(x)*cos(y);
  boundary ({s});

  output_cells (stdout);
  dump (file = "restore.dump", list = {s});
}
