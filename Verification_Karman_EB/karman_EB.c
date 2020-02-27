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
#include "../src_local/output_vtu_foreach.h"

scalar f[];
scalar * tracers = {f};
face vector muv[];

/**
The domain is eight units long, centered vertically. */

int main() {
    L0 = 8.;
    origin (-1.5, -L0/2.);
    N = 512;
    mu = muv;
    TOLERANCE = 1e-8;
    run();
}

/**
We set a constant viscosity corresponding to a Reynolds number of 160,
based on the cylinder diameter (0.125) and the inflow velocity (1). */

event properties (i++)
{
    foreach_face() muv.x[] = fm.x[]*0.125/160.;
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

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

event init (t = 0)
{
    /**
    The domain is the intersection of a channel of width unity and a
    circle of diameter 0.125. */

    vertex scalar phi[];
    foreach_vertex() {
        phi[] = HUGE;
        phi[] = intersection (phi[], sq(x) + sq(y) - sq(0.5));
    }
    boundary ({phi});
    fractions (phi, cs, fs);

    /**
    We set the initial velocity field. */

    foreach()
    u.x[] = cs[];
}

/**
We check the number of iterations of the Poisson and viscous
problems. */

event logfile (i++)
    fprintf (ferr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

/**
We produce animations of the vorticity and tracer fields... */

event movies (i += 4; t <= 15.)
{
scalar omega[], m[];
vorticity (u, omega);
foreach()
        m[] = cs[] - 0.5;
                boundary ({m});
output_ppm (omega, file = "vort.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
        min = -10, max = 10, linear = true, mask = m);
output_ppm (f, file = "f.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
        linear = false, min = 0, max = 1, mask = m);
}

//Output
event vtk_file (t += 0.1){
    char subname[80]; sprintf(subname, "rkEB");
    scalar l[], omega[];
    vorticity (u, omega);
    foreach() {l[] = level;}
    output_vtu_MPI( (scalar *) {cs, f, omega, p, l}, (vector *) {u}, subname, 0 );
}
/**
We adapt according to the error on the embedded geometry, velocity and
tracer fields. */

event adapt (i++) {
    adapt_wavelet ({cs,f,u}, (double[]){1e-2,3e-2,3e-2,3e-2}, 10, 4);
}

/**
## See also

* [Same example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/cylinder.html)
*/