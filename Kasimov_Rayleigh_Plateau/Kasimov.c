
//#include "vtk.h"
//int iteration = 0;
///**
//# Self-similar pinch-off
//
//Here we check whether we can recover the self-similar scalings with
//time of the minimum radius and maximum axial velocity. */
//
//#include "axi.h"
//#include "navier-stokes/centered.h"
//#include "two-phase.h"
//#include "tension.h"
//
//double R = 0.2, epsilon = 0.1;
//
//int main() {
//    rho2 = 0.01;
//    f.sigma = 1.;
//    N = 256;
//    run();
//}
//
//event init (t = 0) {
//    fraction (f, R*(1. + epsilon*cos(pi*x)) - y);
//}
//
//event adapt (i++) {
//    adapt_wavelet ({f}, (double[]){0}, 8); // t < 0.7 ? 8 : 9);
//}
//
//event gfsview (i += 10; t <= 1) {
//static FILE * fp = popen ("gfsview2D -s plateau3.gfv", "w");
//output_gfs (fp);
//}
//
//event logfile (i++) {
//    scalar hy[];
//    foreach() {
//        if (f[] > 1e-3 && f[] < 1 - 1e-3) {
//            coord m = mycs (point, f), fc;
//            double alpha = plane_alpha (f[], m);
//            plane_area_center (m, alpha, &fc);
//            hy[] = y + Delta*fc.y;
//        }
//        else
//            hy[] = nodata;
//    }
//    stats s = statsf (hy);
//    fprintf (stderr, "%g %g %g %g\n", t, s.min, s.max, statsf(u.x).max);
//}


//
//#if _MPI != 1
//event vtk_file (i += 1)
//{
//    char name[80];
//    sprintf(name, "list.vtk.%d", iteration++);
//    printf("name %s", name);
//    FILE * fplist = fopen (name, "w");
//    scalar l[];
//    foreach()
//        l[] = level;
//    output_vtk ({f, u.x, u.y, p, l}, 512, fplist, false );
//}
//#endif



#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "vtk.h"
#include "view.h"
/**
We define the radius of the jet, the initial jet length, the Reynolds
number and the surface tension coefficient. */

#define LEVEL 10
#define length 0.1
#define feps 0
//#define ueps 0.01
#define sharpness 4
/**
The default maximum level of refinement is 10 and the error threshold
on velocity is 0.1. */

int MAXLEVEL = 12;

static int iteration=0;
/**
To impose boundary conditions on a disk we use an auxilliary volume
fraction field *f0* which is one inside the cylinder and zero
outside. We then set an oscillating inflow velocity on the
left-hand-side and free outflow on the right-hand-side. */

scalar f0[];
//u.n[left]  = dirichlet(0);
//u.n[left]  = dirichlet(f0[]*(1. + 0.1*sin (10.*2.*pi*t)));
//u.t[left]  = dirichlet(0);
//#if dimension > 2
//u.r[left]  = dirichlet(0);
//#endif
//p[left]    = neumann(0);
//f[left]    = neumann(0.0);

//u.n[right] = neumann(0);
//u.t[right] = neumann(0);
p[right]   = dirichlet(0);
//f[right]   = neumann(0);


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
    init_grid (1 <<LEVEL);
    origin (0, 0, 0);
    size (1.);
    periodic (top);
    periodic (right);
    /**
    We set the density and viscosity of each phase as well as the
    surface tension coefficient and start the simulation. */

    rho1 = 998.; rho2 = 1.2;
    mu1 = 1.003e-3; mu2 = 1.8e-5;
    f.sigma = 1;
    //G.y = - 1./sq(FROUDE);
    run();
}

/**
## Initial conditions */

double front (double x, double y, double z)
{
//    double phi = HUGE;
    return -pow(2*x-1, sharpness) - pow(60*y-30, sharpness) + pow(0.8,sharpness);
//    return min(phi, 0.2*0.5*(1.0 - tanh((fabs(y-L0/2.0) - 0.02)/0.01)) +0.0-x);

    //return min(phi, -fabs(y-0.5)+(length - x));
    //return min(phi, -0.005*sin(10.0*2.0*pi*y)+(length - x));
}


event init (t = 0) {

    if (!restore (file = "restart")) {

        /**
        We use a static refinement down to *maxlevel* in a cylinder 1.2
        times longer than the initial jet and twice the radius. */
        int it = 0;
        do {
            it++;
            fraction(f0, front (x, y, z));
            if (it>=10) printf("WARNING: does not converge... ");
        }while (adapt_wavelet({f0}, (double []){feps},
                              maxlevel = MAXLEVEL, minlevel = 5).nf != 0 && it <= 10);


        /**
        We then use this to define the initial jet and its velocity. */

        foreach() {
            f[] = f0[];
            u.x[] = 0;
        }
        boundary ({f,u.x});
    }
}

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

//event movie (t += 1e-1; t <= 10)
//{
//#if dimension == 2
//  scalar omega[];
//    vorticity (u, omega);
//  view (tx = 0.5);
//  clear();
//  draw_vof ("f");
//  squares ("omega", linear = true, spread = 10);
//  box ();
//#else // 3D
//    scalar pid[];
//    foreach()
//    pid[] = fmod(pid()*(npe() + 37), npe());
//    boundary ({pid}); // not used for the moment
//    view (camera = "iso",
//          fov = 14.5, tx = -0.418, ty = 0.288,
//          width = 1600, height = 1200);
//    clear();
//    draw_vof ("f");
//#endif // 3D
//    save ("movie.mp4");
//}

//Output
#include "output_fields/output_vtu_foreach.h"
event vtk_file (t += 0.01; t <= 10)
{
    int nf = iteration;
    scalar l[];
    foreach()
    l[] = level;

    char name[80], subname[80];
    FILE *fp;
    sprintf(name, "hs_%4.4d_n%3.3d.vtu", nf, pid());
    fp = fopen(name, "w");

    output_vtu_bin_foreach((scalar *) {f, rhov, p, l}, (vector *) {u}, 64, fp, false);
    fclose(fp);
    @if _MPI
    if (pid() == 0) {
        sprintf(name, "hs_%4.4d.pvtu", nf);
        sprintf(subname, "hs_%4.4d", nf);
        fp = fopen(name, "w");
        output_pvtu_bin((scalar *) {f, rhov, p, l}, (vector *) {u}, 64, fp, subname);
        fclose(fp);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    @endif

    iteration++;
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



/**
## Mesh adaptation

We adapt the mesh according to the error on the volume fraction field
and the velocity. */

event adapt (i++) {
    adapt_wavelet ({f}, (double[]){feps}, MAXLEVEL);
}

