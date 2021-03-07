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

#define FILTERED
#define RELATIVE_RESIDUAL
#define EPS_MAXA 2
#define TURN_ON_TRACER 0
#include "embed.h"
#include "navier-stokes/centered.h"
// #include "navier-stokes/perfs.h"
#include "../src_local/output_vtu_foreach.h"
#include "../src_local/utils-weugene.h"
#if TURN_ON_TRACER == 1
    #include "tracer.h"
    scalar f[];
    scalar * tracers = {f};
#endif
face vector muv[];
scalar omega[], divu[];
int maxlevel = 11;
int minlevel = 4;
double xx0 = 0, rad = 0.5, Ldomain=50, RE=40.;
/**
The fluid is injected on the left boundary with a unit velocity. The
tracer is injected in the lower-half of the left boundary. An outflow
condition is used on the right boundary. */

u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

#if TURN_ON_TRACER == 1
    f[left]   = dirichlet(y < 0);
#endif
/**
The top and bottom walls are free-slip and the cylinder is no-slip. */

u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

void soild_fs_embed(scalar cs, face vector fs, double t){
    vertex scalar phi[];
    foreach_vertex() {
        phi[] = HUGE;
        phi[] = intersection(phi[], ( sq(x - xx0) + sq(y) - sq(rad)));
    }
    boundary ({phi});
    fractions (phi, cs, fs);
//    fraction (fs, sq(rad) - sq(x - xx0) - sq(y));
}

int main(int argc, char * argv[]) {
    if (argc > 1) {
        maxlevel = atoi(argv[1]); //convert from string to float
    }
    size (Ldomain);
    origin (-0.3*Ldomain, -0.5*Ldomain);
    DT = 1e-8;
    CFL=0.4;
    TOLERANCE = 1e-8;
    RELATIVE_RES_TOLERANCE = 0.1;
    NITERMAX=30;
    N = 512;
    mu = muv;

    run();
}

/**
We set a constant viscosity corresponding to a Reynolds number of 40, 100,
based on the cylinder diameter (1) and the inflow velocity (1). */

event properties (i++)
{
    foreach_face() muv.x[] = 2*rad*fm.x[]/RE;
}

event init (t = 0)
{
    /**
    The domain is the intersection of a channel of width unity and a
    circle of diameter 0.125. */

//    vertex scalar phi[];
//    foreach_vertex() {
//        phi[] = HUGE;
//        phi[] = intersection (phi[], sq(x) + sq(y) - sq(0.5));
//    }
//    boundary ({phi});
//    fractions (phi, cs, fs);
    if (!restore (file = "restart")) {
        int it = 0;
#if TURN_ON_TRACER == 1
        do {
            soild_fs_embed (cs, fs, 0);
            foreach() f[] = 0;
            boundary({f, cs});
        }while (adapt_wavelet({cs, f}, (double []){1e-5, 1e-5}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
#else
        do {
            soild_fs_embed (cs, fs, 0);
            boundary({cs});
        }while (adapt_wavelet({cs}, (double []){1e-5}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
#endif
        /**
        We set the initial velocity field. */
        foreach() u.x[] = cs[];
        boundary(all);
        event("vtK_file");
    }


}

event set_dtmax (i++) {
    if (i<=100) {
        NITERMIN=100;
        NITERMAX=150;
    }else{
        NITERMIN=10;
        NITERMAX=30;
    }
    DT *= 1.05;
    DT = min(DT, CFL*Ldomain/pow(2, maxlevel+3));
    fprintf(ferr, "set_dtmax: tnext= %g Dt= %g", tnext, DT);
}

/**
We check the number of iterations of the Poisson and viscous
problems. */

event logfile (i++){
    foreach() {
        divu[] = 0;
        foreach_dimension() divu[] += (uf.x[1]-uf.x[])/Delta;
    }
    double Linfu = -10;
    foreach( reduction(max:Linfu) ){
        if (fabs(divu[]) > Linfu) Linfu = fabs(divu[]);
    }
    fprintf (ferr, "i=%d t=%g dt=%g iter_p=%d iter_u=%d div u=%g \n", i, t, dt, mgp.i, mgu.i, Linfu);
}
/**
We produce animations of the vorticity and tracer fields... */

//event movies (i += 4; t <= 15.)
//{
//    scalar omega[], m[];
//    vorticity (u, omega);
//    foreach()
//            m[] = cs[] - 0.5;
//    boundary ({m});
//    output_ppm (omega, file = "vort.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
//        min = -10, max = 10, linear = true, mask = m);
//#if TURN_ON_TRACER == 1
//    output_ppm (f, file = "f.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
//        linear = false, min = 0, max = 1, mask = m);
//#endif
//}

//Output
//event vtk_file (i++){
event vtk_file (t += 1){
    char subname[80]; sprintf(subname, "rkEB");
    scalar l[], omega[];
    vorticity (u, omega);
    foreach() {l[] = level;}
#if TURN_ON_TRACER == 1
    output_vtu_MPI( (scalar *) {cs, f, omega, p, l}, (vector *) {u}, subname, 0 );
#else
    output_vtu_MPI( (scalar *) {cs,    omega, p, l}, (vector *) {u}, subname, 0 );
#endif
}
/**
We adapt according to the error on the embedded geometry, velocity and
tracer fields. */
#if TURN_ON_TRACER == 1
#define ADAPT_SCALARS {cs, f, u}
#define ADAPT_EPS_SCALARS {1e-5, 1e-5, 1e-2, 1e-2}
#else
#define ADAPT_SCALARS {cs, u}
#define ADAPT_EPS_SCALARS {1e-5, 1e-3, 1e-3}
#endif
event adapt (i++){
    double eps_arr[] = ADAPT_EPS_SCALARS;
    MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
}

event stop (t = 500);
/**
## See also

* [Same example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/cylinder.html)
*/
