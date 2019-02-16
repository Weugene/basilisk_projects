/**
# Simple test of Basilisk View

Ran in parallel on four MPI processes. */

#include "fractions.h"
#include "view.h"

int main() {

  /**
  We first define a volume fraction field. */
  
  init_grid (16);
  origin (-0.5,-0.5,-0.5);
  scalar f[];
  fraction (f, sq(x) + sq(y) + sq(z) - sq(0.3));

  /**
  Then display it using Basilisk view functions. 
  
  <table>
  <tr>
  <td>![2D Basilisk view](view/out.png)</td>
  <td>![3D Basilisk view](view.3D/out.png)</td>
  </tr>
  <tr>
  <td>2D Basilisk view</td>
  <td>3D Basilisk view</td>
  </tr>
  </table>
  */
  
  view (width = 400, height = 400);
  box();
  draw_vof ("f");
  cells();
  squares ("f", min = 0, max = 1);
  save ("out.png");

  /**
  A few more files just for testing. */
  
  output_facets (f, stderr);
  dump (file = "dump");
}
