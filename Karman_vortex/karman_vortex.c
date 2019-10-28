/**
# Bénard–von Kármán Vortex Street for flow around a cylinder at Re=160

An example of 2D viscous flow around a simple solid boundary. Fluid is
injected to the left of a channel bounded by solid walls with a slip
boundary condition. A passive tracer is injected in the bottom half of
the inlet.

![Animation of the vorticity field.](karman/vort.mp4)(loop)

![Animation of the tracer field.](karman/f.mp4)(loop)

We use the centered Navier-Stokes solver, with embedded boundaries and
advect the passive tracer *f*. */

#include "embed.h"
#include "navier-stokes/centered.h"
// #include "navier-stokes/perfs.h"
#include "tracer.h"

scalar f[];
scalar * tracers = {f};
face vector muv[];
#define MAXLEVEL 12
/**
The domain is eight units long, centered vertically. */

int main() {
  L0 = 8.;
  origin (-0.5, -L0/2.);
  N = 1024;
  mu = muv;
  run(); 
}

/**
We set a constant viscosity corresponding to a Reynolds number of 160,
based on the cylinder diameter (0.125) and the inflow velocity (1). */

event properties (i++)
{
  foreach_face()
    muv.x[] = fm.x[]*0.125/160.;
}

/**
The fluid is injected on the left boundary with a unit velocity. The
tracer is injected in the lower-half of the left boundary. An outflow
condition is used on the right boundary. */

u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(y < 0);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
The top and bottom walls are free-slip and the cylinder is no-slip. */

u.n[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);
u.t[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);

event init (t = 0)
{
//  refine (sq(x) + sq(y) < sq(0.125) && level < MAXLEVEL);
  /**
  The domain is the intersection of a channel of width unity and a
  circle of diameter 0.125. */

  vertex scalar phi[];
  foreach_vertex() {
    phi[] = intersection (0.5 - y, 0.5 + y);
    phi[] = intersection (phi[], sq(x) + sq(y) - sq(0.125/2.));
  }
  boundary ({phi});
  fractions (phi, cs, fs);

  /**
  We set the initial velocity field. */

  foreach(){
    u.x[] = cs[] ? 1. : 0.;
    f[] = 0;
  }
  boundary({f, u});
  event("end_timestep");
}



/**
We check the number of iterations of the Poisson and viscous
problems. */

event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
scalar omega[];
//Output
#include "../src_local/output_vtu_foreach.h"
event end_timestep (i++){
    char subname[80]; sprintf(subname, "kk");
    scalar l[];
    vorticity (u, omega);
    foreach() l[] = level;
    output_vtu_MPI( (scalar *) {l, omega, p, cs, f}, (vector *) {u}, subname);
}

/**
We adapt according to the error on the embedded geometry, velocity and
tracer fields. */

#define ADAPT_SCALARS {cs, omega}
#define ADAPT_EPS_SCALARS {1e-3,1e-3,1e-3,1e-3}
event adapt (i++) {
  fprintf(stderr, "before");
  double eps_arr[] = ADAPT_EPS_SCALARS;
//  MinMaxValues(ADAPT_SCALARS, eps_arr);
  adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = 10, minlevel = 4);
  fprintf(stderr, "after");
}



event stop(t = 10);
/**
## See also

* [Same example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/cylinder.html)
*/
/**
We produce animations of the vorticity and tracer fields... */

//event movies (i += 4; t <= 15.)
//{
//  scalar omega[], m[];
//  vorticity (u, omega);
//  foreach()
//    m[] = cs[] - 0.5;
//  boundary ({m});
//  output_ppm (omega, file = "vort.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
//	      min = -10, max = 10, linear = true, mask = m);
//  output_ppm (f, file = "f.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
//	      linear = false, min = 0, max = 1, mask = m);
//}