/**
# Self-similar pinch-off

Here we check whether we can recover the self-similar scalings with
time of the minimum radius and maximum axial velocity. */

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

#define MAXLEVEL 10
double R = 0.05, epsilon = 0.1;
int iteration = 0, itermain = 0;
int main(int argc, char * argv[]) {
    if (argc > 1)
        itermain = atoi (argv[1]);
    periodic(right);
    rho2 = 0.01;
    f.sigma = 1.;
    N = 1 << MAXLEVEL;
    run();
}

event init (t = 0) {
    if (!restore (file = "restart")) {
//    -pow(2*x-1, sharpness) - pow(60*y-30, sharpness) + pow(0.8,sharpness)
        fraction(f, R * (1. + epsilon * cos(11*pi * x )) - y);
    }else{ iteration = itermain;}
}

event adapt (i++) {
    adapt_wavelet ({f}, (double[]){0}, MAXLEVEL); // t < 0.7 ? 8 : 9);
}

/*event gfsview (i += 10; t <= 1) {
static FILE * fp = popen ("gfsview2D -s plateau3.gfv", "w");
output_gfs (fp);
}*/

event logfile (i++) {
    scalar hy[];
    foreach() {
        if (f[] > 1e-3 && f[] < 1 - 1e-3) {
            coord m = mycs (point, f), fc;
            double alpha = plane_alpha (f[], m);
            plane_area_center (m, alpha, &fc);
            hy[] = y + Delta*fc.y;
        }
        else
            hy[] = nodata;
    }
    stats s = statsf (hy);
    fprintf (stderr, "%g %g %g %g\n", t, s.min, s.max, statsf(u.x).max);
}


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

event snapshot (i += 10000) {
    char name[80];
    sprintf (name, "dump-%d", i);
    scalar pid[];
    foreach()
        pid[] = fmod(pid()*(npe() + 37), npe());
    boundary ({pid});
    dump (name);
}




