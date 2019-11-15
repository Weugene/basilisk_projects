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

//#include "embed.h"
//#include "navier-stokes/centered.h"
vector tmp[];
scalar fs[];    //added
scalar f[];
#define BRINKMAN_PENALIZATION 1
#define DEBUG_BRINKMAN_PENALIZATION
#include "../src_local/centered-weugene.h"

#include "fractions.h" //added
// #include "navier-stokes/perfs.h"
#include "tracer.h"

#define Re 160.
#define diam 0.125
#define Lch 0.5
#define MAXLEVEL 11

scalar * tracers = {f};
face vector muv[];
bool flag = true;
/**
The domain is eight units long, centered vertically. */

int main() {
    L0 = 8.;
    origin (-0.5, -L0/2.);
    N = 1024;
    mu = muv;
    eta_s = 1e-6;
    TOLERANCE = 1e-5;
//    nu_s = sq(L0/pow(2.,MAXLEVEL))/eta_s/100;
//    stokes = true;
    run();
}

/**
We set a constant viscosity corresponding to a Reynolds number of 160,
based on the cylinder diameter (0.125) and the inflow velocity (1). */

event properties (i++)
{
    foreach_face()
    muv.x[] = fm.x[]*diam/Re;
//    foreach() foreach_dimension() target_U.x[] = 0;
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

//u.n[embed] = fabs(y) > Lch/2 ? neumann(0.) : dirichlet(0.);
//u.t[embed] = fabs(y) > Lch/2 ? neumann(0.) : dirichlet(0.);

event init (t = 0)
{

    /**
    The domain is the intersection of a channel of width unity and a
    circle of diameter diam. */

//    vertex scalar phi[];
//    face vector face_fs[];
//    foreach_vertex() {
////        phi[] = intersection (L0 - y, 0.5 + y);
//        phi[] = (sq(x) + sq(y) - sq(diam/2.) <= 0) ? 1 : -1;
//    }
//    boundary ({phi});
//    fractions (phi, fs, face_fs);
    foreach() fs[] = (sq(x) + sq(y) - sq(diam/2.) <= 0) ? 1 : 0;

    /**
    We set the initial velocity field. */

    foreach()
    u.x[] = fs[] ? 0. : 1.;
    event("properties");
    event("end_timestep");
}

event velocity_correction(i++){
    foreach() {
        fs[] = (sq(x) + sq(y) - sq(diam/2.) <= 0) ? 1 : 0;
        foreach_dimension() {u.x[] *= 1 - fs[];}
    }
    foreach_face() {
        uf.x[] *= 1 - fs[];
    }
    if (flag) event ("end_timestep");
    flag = false;
}
/**
We check the number of iterations of the Poisson and viscous
problems. */

event logfile (i++)
fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

/**
We produce animations of the vorticity and tracer fields... */

event movies (t += 0.1; t <= 15.)
{
scalar omega[], m[];
vorticity (u, omega);
foreach()
m[] = fs[] - 0.5;
boundary ({m});
output_ppm (omega, file = "vort.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
        min = -10, max = 10, linear = true, mask = m);
output_ppm (f, file = "f.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
        linear = false, min = 0, max = 1, mask = m);
}


//Output
#include "../src_local/output_vtu_foreach.h"
event end_timestep (t += 0.1){
    char subname[80]; sprintf(subname, "br");
    scalar l[], omega[];
//    vector ve[]; scalar pe[]; exact_solution(ve, pe);
    vorticity (u, omega);
    foreach() l[] = level;
//    output_vtu_MPI( (scalar *) {l, omega, p, f}, (vector *) {u}, subname, 0);
    output_vtu_MPI( (scalar *) {l, omega, p, pf, f, fs}, (vector *) {u, uf, mu, total_rhs, dbp, tmp}, subname, 0);
    flag = true;
}

/**
We adapt according to the error on the embedded geometry, velocity and
tracer fields. */

event adapt (i++) {
    adapt_wavelet ({fs,u,f}, (double[]){1e-2,3e-2,3e-2,3e-2}, MAXLEVEL, 4);
}

/**
## See also

* [Same example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/cylinder.html)
*/