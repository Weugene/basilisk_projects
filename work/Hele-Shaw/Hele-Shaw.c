/**
# Atomisation of a pulsed liquid jet

A dense cylindrical liquid jet is injected into a stagnant lighter
phase (density ratio 1/27.84). The inflow velocity is modulated
sinusoidally to promote the growth of primary shear
instabilities. Surface tension is included and ultimately controls the
characteristic scale of the smallest droplets.

We solve the two-phase Navier--Stokes equations with surface
tension. We need the *tag()* function to count the number of
droplets. We generate animations online using Basilisk View. */
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "embed.h"

#include "tag.h"
#include "view.h"
//ONLY for 3D case. for 2D case it doesn't work
#include "fractions_output.h"
#include "output_fields/output_vtu_foreach.h"
/**
We define the radius of the jet, the initial jet length, the Reynolds
number and the surface tension coefficient. */

#define LEVEL 6
#define length 0.1
#define height 0.05
#define rad 0.05
#define feps 0.01
#define ueps 0.01
#define sharpness 4
/**
The default maximum level of refinement is 10 and the error threshold
on velocity is 0.1. */

int MAXLEVEL = 8;
static int iteration=0;

/**
To impose boundary conditions on a disk we use an auxilliary volume
fraction field *f0* which is one inside the cylinder and zero
outside. We then set an oscillating inflow velocity on the
left-hand-side and free outflow on the right-hand-side. */

scalar f0[], hole[];
scalar channel[];
//face vector muv[];
u.n[left]  = neumann(0);
u.t[left]  = neumann(0);
u.r[left]  = neumann(0);
p[left]    = neumann(0);
f[left]    = neumann(0);

u.n[right] = neumann(0);
u.t[right] = neumann(0);
u.r[right] = neumann(0);
p[right]   = dirichlet(0);
f[right]   = neumann(0);

u.n[bottom]  = neumann(0);
u.t[bottom]  = neumann(0);
u.r[bottom]  = neumann(0);
p[bottom]    = neumann(0);
f[bottom]    = neumann(0);

u.n[top] = neumann(0);
u.t[top] = neumann(0);
u.r[top] = neumann(0);
p[top]   = neumann(0);
f[top]   = neumann(0);
//front does not work
u.n[front]  = dirichlet(0);
u.t[front]  = dirichlet(0);
u.r[front]  = dirichlet(0);
p[front]    = neumann(0);
f[front]    = neumann(0);

//u.n[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);
//u.t[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);
u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);
u.r[embed] = dirichlet (0.);
//p[embed] = neumann (0.);
//f[embed] = dirichlet (0.);


u.n[back]  = dirichlet(hole[]);
u.t[back]  = dirichlet(0);
u.r[back]  = dirichlet(0);
p[back]    = neumann(0);
f[back]    = dirichlet(hole[]);


/**
The program can take two optional command-line arguments: the maximum
level and the error threshold on velocity. */

int main (int argc, char * argv[])
{
    if (argc > 1)
        MAXLEVEL = atoi (argv[1]);
//    if (argc > 2)
//        ueps = atof (argv[2]);

    /**
    The initial domain is discretised with $64^3$ grid points. We set
    the origin and domain size. */
    init_grid (1 << LEVEL);
    origin (0, 0, 0);
    size (1.);
    stokes = true;
    /**
    We set the density and viscosity of each phase as well as the
    surface tension coefficient and start the simulation. */

    rho1 = 1000; rho2 = 100;
    mu1 = 1.e-3; mu2 = 2.0e-5;
    f.sigma = 1;
    //mu = muv;
    //G.y = - 1./sq(FROUDE);
    run();
}

//event properties (i++)
//{
//    foreach_face()
//    muv.x[] = fm.x[]/300.;
//}
/**
## Initial conditions */

double front_geom (double x, double y, double z)
{
    return rad*rad - pow(x - 0.5*L0, 2) - pow(y - 0.5*L0, 2);
//    double phi = HUGE;
    //return -pow(2*x-1, sharpness) - pow(20*y-10, sharpness) + pow(0.8,sharpness);
//    return min(phi, 0.2*0.5*(1.0 - tanh((fabs(y-L0/2.0) - 0.02)/0.01)) +0.0-x);

    //return min(phi, -fabs(y-0.5)+(length - x));
    //return min(phi, -0.005*sin(10.0*2.0*pi*y)+(length - x));
}

double channel_geom (double x, double y, double z)
{
    return height - z;
}

double hole_geom (double x, double y, double z)
{
    double phi = rad*rad - pow(x - 0.5*L0, 2) - pow(y - 0.5*L0, 2);
    return min(phi, height - z);
}

event init (t = 0) {
    if (!restore (file = "restart")) {
        /**
        We use a static refinement down to *maxlevel* in a cylinder 1.2
        times longer than the initial jet and twice the radius. */
        int it = 0;

        do {
            it++;
            printf("channel");
            fraction(channel, channel_geom (x, y, z));
            printf("hole");
            fraction(hole, hole_geom (x, y, z));
            if (it>=10) printf("WARNING: does not converge... ");
        }while (adapt_wavelet({channel, hole}, (double []){feps, feps},
                              maxlevel = MAXLEVEL, minlevel = 3).nf != 0 && it <= 10);

        /**
        We then use this to define the initial jet and its velocity. */

        foreach() {
            f[] = hole[];
            f0[]= hole[];
            foreach_dimension()
                u.x[] = 0;
        }
        vertex scalar phi[];
        foreach_vertex()
            phi[] = height - z - 0.001;
        boundary ({phi});
        fractions (phi, cs, fs);
        fractions (channel, cs, fs);

        u.n[embed] = dirichlet (0.);
        u.t[embed] = dirichlet (0.);
        u.r[embed] = dirichlet (0.);
        f[embed]   = dirichlet (0.);

        u.n[back] = dirichlet (hole[]*min(100.0*t, 1.0));
        u.t[back] = dirichlet (0.);
        u.r[back] = dirichlet (0.);
        f[back]   = dirichlet (hole[]);

        boundary ({f, u.x, u.y, u.z, p});

        int nf = 0;
        if (iteration>0) {
            char name[80], subname[80];
            FILE *fp;
            sprintf(name, "hs_%4.4d_n%3.3d.vtu", nf, pid());
            fp = fopen(name, "w");
            scalar l[];
            foreach()
            l[] = level;
            output_vtu_bin_foreach((scalar *) {f, hole, channel, l, p, rhov, mu}, (vector *) {u, mu, fm}, 64, fp, false);
            fclose(fp);
            @if _MPI
            if (pid() == 0) {
                sprintf(name, "hs_%4.4d.pvtu", nf);
                sprintf(subname, "hs_%4.4d", nf);
                fp = fopen(name, "w");
                output_pvtu_bin((scalar *) {f, hole, channel, l, p, rhov, mu}, (vector *) {u, mu, fm}, 64, fp, subname);
                fclose(fp);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            @endif
        }
    }
//    char name[80];
//    sprintf(name, "list.ply.%d", iteration++);
//    printf("name %s", name);
//    FILE * fp = fopen (name, "w");
//    output_ply(f, fp);
}

//event velocity (i++) {
//    foreach()
//    foreach_dimension()
//    u.x[] = channel[]*u.x[];//velocity = 0 in the solid
//    boundary ((scalar *){u});
//}

/**
## Outputs

We log some statistics on the solver. */

event logfile (i++) {
    if (i == 0)
        fprintf (ferr,
                 "t dt mgp.i mgpf.i mgu.i grid->tn perf.t perf.speed\n");
    fprintf (ferr, "%g %g %d %d %d %ld %g %g\n",
             t, dt, mgp.i, mgpf.i, mgu.i,
             grid->tn, perf.t, perf.speed);
}

/**
We generate an animation using Basilisk View. */

event movie (t += 1e-1; t <= 10)
{
#if dimension == 2
  scalar omega[];
    vorticity (u, omega);
  view (tx = 0.5);
  clear();
  draw_vof ("f");
  squares ("omega", linear = true, spread = 10);
  box ();
#else // 3D
    scalar pid[];
    foreach()
    pid[] = fmod(pid()*(npe() + 37), npe());
    boundary ({pid}); // not used for the moment
    view (camera = "iso", fov = 14.5, tx = 0.5, ty = 0.5);
    clear();
    draw_vof ("f");
#endif // 3D
    save ("movie.mp4");
}

event vtk_file (i += 1)
{
    int nf = ++iteration;
    if (iteration>0) {
        char name[80], subname[80];
        FILE *fp;
        sprintf(name, "hs_%4.4d_n%3.3d.vtu", nf, pid());
        fp = fopen(name, "w");
        scalar l[];
        foreach()
            l[] = level;
        output_vtu_bin_foreach((scalar *) {f, hole, channel, l, p, rhov}, (vector *) {u, mu, fm}, 64, fp, false);
        fclose(fp);
        @if _MPI
        if (pid() == 0) {
            sprintf(name, "hs_%4.4d.pvtu", nf);
            sprintf(subname, "hs_%4.4d", nf);
            fp = fopen(name, "w");
            output_pvtu_bin((scalar *) {f, hole, channel, l, p, rhov}, (vector *) {u, mu, fm}, 64, fp, subname);
            fclose(fp);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        @endif
    }
}

/**
We save snapshots of the simulation at regular intervals to
restart or to post-process with [bview](/src/bview). */

event snapshot (i += 1000) {
    char name[80];
    sprintf (name, "dump-%d", i);
    scalar pid[];
    foreach()
        pid[] = fmod(pid()*(npe() + 37), npe());
    boundary ({pid});
    dump (name);
}

//event fraction_output (i += 1)
//{
//    char name[80];
//    sprintf(name, "list.ply.%d", iteration++);
//    printf("name = %s", name);
//    FILE * fplist = fopen (name, "w");
//    output_ply(f, fplist);
//}



/**
## Mesh adaptation

We adapt the mesh according to the error on the volume fraction field
and the velocity. */

event adapt (i++) {
    adapt_wavelet ({f, channel, u.x, u.y, u.z}, (double[]){feps, feps, ueps, ueps, ueps}, MAXLEVEL);
}
