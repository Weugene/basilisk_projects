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
//#include "navier-stokes/centered.h"
scalar f[];
#include "../src_local/centered-weugene.h"
#include "fractions.h" //added
// #include "navier-stokes/perfs.h"
#include "tracer.h"

#define Re 1.
#define diam 0.125
#define rad (0.5*diam)
#define Lch 0.5
#define MAXLEVEL 11

scalar * tracers = {f};
face vector muv[];
bool flag = true;
/**
The domain is eight units long, centered vertically. */
void exact_solution(vector ve, scalar pe){
  double theta, r, vr, vth;
  foreach() {
    r = sqrt(sq(x) + sq(y));
    theta = atan2(y, x + SEPS);
#if 1
    vr  =   (r > rad) * (1 - sq(rad/(r + SEPS))) * cos(theta);
    vth = - (r > rad) * (1 + sq(rad/(r + SEPS))) * sin(theta);
    ve.x[] = vr * cos(theta) - vth * sin(theta);
    ve.y[] = vr * sin(theta) + vth * cos(theta);
    pe[] = 0.5 * rho[] * (1 - sq(ve.x[]) - sq(ve.y[]));
#endif
  }
}



int main() {
    L0 = 8.;
    origin (-0.5, -L0/2.);
    N = 64;
    mu = muv;
    TOLERANCE = 1e-5;
    DT=1e-4;
//    nu_s = sq(L0/pow(2.,MAXLEVEL))/eta_s/100;
    stokes = true;
    run();
}
/**
The fluid is injected on the left boundary with a unit velocity. The
tracer is injected in the lower-half of the left boundary. An outflow
condition is used on the right boundary. */

u.n[left]  = dirichlet(1.);
//uf.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
//pf[left]   = neumann(0.);
f[left]    = dirichlet(y < 0);

u.n[right] = neumann(0.);
//uf.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
//pf[right]  = dirichlet(0.);

/**
The top and bottom walls are free-slip and the cylinder is no-slip. */

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
/**
We set a constant viscosity corresponding to a Reynolds number of 160,
based on the cylinder diameter (0.125) and the inflow velocity (1). */

event properties (i++)
{
    foreach_face() muv.x[] = fm.x[]*diam/Re;
//    foreach() foreach_dimension() target_U.x[] = 0;
}


event init (t = 0)
{
    refine (((sq(x) + sq(y) < sq(1.5*rad)) || ((x< -0.3) && (y < rad))) && (y < 2*rad) && (level < MAXLEVEL));

    vertex scalar phi[];
    foreach_vertex() {
        phi[] = (sq(x) + sq(y) < sq(diam/2.) ) ? -1 : 1;
    }
    boundary ({phi});
    fractions (phi, cs, fs);
    foreach() f[] = (x<-0.4)*(y<0);

    foreach() u.x[] = 1.0 - cs[];
    foreach_face() uf.x[] = 1.0 - fs.x[];
    event("properties");
    event("end_timestep");
}

event velocity_correction(i++){
//    foreach() {
//        cs[] = (sq(x) + sq(y) - sq(diam/2.) <= 0) ? 1 : 0;
//        foreach_dimension() {u.x[] *= 1 - cs[];}
//    }
//    foreach_face() {
//        uf.x[] *= 1 - fs.x[];
//    }
    if (flag) event ("end_timestep");
    flag = false;
}
/**
We check the number of iterations of the Poisson and viscous
problems. */

event logfile (i++)
fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

//Output
#include "../src_local/output_vtu_foreach.h"
event end_timestep (t += 1e-1; t <= 15.){
//event end_timestep (i += 10; t <= 15.){
//    if (i<=1) return 0;
    char subname[80]; sprintf(subname, "br");
    scalar l[], omega[];
//    vector ve[]; scalar pe[]; exact_solution(ve, pe);
    vorticity (u, omega);
    foreach() l[] = level;
//    output_vtu_MPI( (scalar *) {l, omega, p, f}, (vector *) {u}, subname, 0);
    output_vtu_MPI( (scalar *) {l, omega, p, pf, f, cs}, (vector *) {u, uf}, subname, 0);
    flag = true;
}

/**
We adapt according to the error on the embedded geometry, velocity and
tracer fields. */

event adapt (i++) {
    adapt_wavelet ({cs,u,f}, (double[]){1e-2,3e-2,3e-2,3e-2}, MAXLEVEL, 4);
}
