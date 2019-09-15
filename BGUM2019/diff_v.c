/**
# Tutorial part 2: Variable diffusivity and an adaptive grid

Here we extend the code from the [more simple diffusion
example](diff.c) and set a variable diffusivity and use an adaptive
grid. Please note the note in the corresponding event.

Results:

![The evolution of the scalar field `s`](diff_v/s.mp4)

![The level of refinement](diff_v/level.mp4)
*/
#include "diffusion.h"
#include "run.h" 

scalar s[];

int main() {
  DT = 0.05;
  run();
}

event init (t = 0) {
  foreach() 
    s[] = exp(-(sq(x - 0.5) + sq(y - 0.5))*10.); 
}

event mov (t += 0.1) {
  output_ppm (s, file = "s.mp4", n = 256, min = -1, max = 1);
  scalar lev[];
  foreach() 
    lev[] = level; //`level` is available in the grid iterator
  output_ppm (lev, file = "level.mp4", n = 256, max = 6);
}

event diff (i++) {
  dt = dtnext (DT);
  /**
     We define a variable diffusivity field (`kap`) and call the
     `boudary()` function in order to ensure its proper difiniton near
     resolutions boundaries and on the various lower levels. [see
      this documentation page](/sandbox/Antoonvh/The_Tree_Grid_Structure_in_Basilsik). Note
     for those who attended the tutorial: I originally *forgot* to
     include this important step...
  */
  face vector kap[];
  foreach_face(x)
    kap.x[] = x/100.;
  foreach_face(y)
    kap.y[] = 0.;
  boundary ((scalar*) {kap});
  diffusion (s, dt, kap);
}
/**
   We "switch on" adaptivity with this event:
*/
event adapt (i++) 
  adapt_wavelet ({s}, (double[]){0.01}, 6);

event stop (t = 10);