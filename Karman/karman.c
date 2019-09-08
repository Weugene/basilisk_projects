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
#include "navier-stokes/centered.h"
// #include "navier-stokes/perfs.h"
#include "tracer.h"

scalar f[];
scalar * tracers = {f};
face vector muv[];
face vector av[];
scalar cs[];//an analogue is
/**
The domain is eight units long, centered vertically. */

int main() {
    L0 = 8.;
    origin (-0.5, -L0/2.);
    N = 2048;
    mu = muv;
    a = av;
    run();
}

/**
We set a constant viscosity corresponding to a Reynolds number of 160,
based on the cylinder diameter (0.125) and the inflow velocity (1). */

event properties (i++)
{
    foreach_face()
        muv.x[] = fm.x[]*0.125/160.*(1+(x>2&& x<3.5)*0);
//        muv.x[] = fm.x[]*0.125/160.*(1+(x>2&& x<5)*100);
}

/**
The fluid is injected on the left boundary with a unit velocity. The
tracer is injected in the lower-half of the left boundary. An outflow
condition is used on the right boundary. */

u.n[left]  = dirichlet(0.1*t);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(y < 0);//a lower half

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
The top and bottom walls are free-slip and the cylinder is no-slip. */

//u.n[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);
//u.t[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);
//vertex scalar phi[];
event init (t = 0)
{
    /**
    The domain is the intersection of a channel of width unity and a
    circle of diameter 0.125. */
//    foreach_vertex() {
//        phi[] = intersection (0.5 - y, 0.5 + y);
//        phi[] = intersection (phi[], sq(x) + sq(y) - sq(0.125/2.));
//    }
//    boundary ({phi});
//    fractions (phi, cs, fs);
    foreach(){
        cs[] =(sq(x) + sq(y) - sq(0.125/2.)>0)?1:0;
    }
    /**
    We set the initial velocity field. */

    foreach()
    {
        u.x[] = 0.;
        f[] = 0;
    }
}

/**
We check the number of iterations of the Poisson and viscous
problems. */

event logfile (i++)
fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

event acceleration (i++) {
    foreach_face()
        av.x[] = (cs[]/(0.01))*((u.x[] + u.x[-1])/2. );
//    foreach_face()
//        av.x[] -= cs[]*(u.x[]);
//    boundary ((scalar *){av});
}
event stability (i++,last) {
    dt = dtnext (stokes ? dtmax : timestep (fabs(uf)+fabs(dt*av.x), dtmax));
}
/**
We produce animations of the vorticity and tracer fields... */

event movies (t+=0.01; i <= 5000.)
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
static int iteration=0;
#include "output_fields/output_vtu_foreach.h"
event vtk_file (t += 0.05)
{
    int nf = iteration;
    scalar l[];
    foreach()
        l[] = level;

    char name[80], subname[80];
    FILE *fp;
    sprintf(name, "hrhs_%4.4d_n%3.3d.vtu", nf, pid());
    fp = fopen(name, "w");

    output_vtu_bin_foreach((scalar *) {l, f, cs}, (vector *) {u, mu}, 64, fp, false);
    fclose(fp);
    @if _MPI
        if (pid() == 0) {
            sprintf(name, "hrhs_%4.4d.pvtu", nf);
            sprintf(subname, "hrhs_%4.4d", nf);
            fp = fopen(name, "w");
            output_pvtu_bin((scalar *) {l, f, cs}, (vector *) {u, mu}, 64, fp, subname);
            fclose(fp);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    @endif
    fprintf (ferr, "iteration: %d\n", iteration); fflush (ferr);
    iteration++;
}




#if DUMP
event snapshot (i += 1000)
//event snapshot (i++)
{
  char name[80];
  sprintf (name, "dump-%d", i);
//  scalar l2[];

//  lambda2 (u, l2);
  dump (file = name);
}
#endif

/**
We adapt according to the error on the embedded geometry, velocity and
tracer fields. */

event adapt (i++) {
    adapt_wavelet ({cs,u,f}, (double[]){1e-3,1e-3,1e-3,1e-3}, 11, 4);
}

/**
## See also

* [Same example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/cylinder.html)
*/