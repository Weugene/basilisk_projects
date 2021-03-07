#include "navier-stokes/centered.h"
#include "navier-stokes/double-projection.h"
#define RAD (pow(pow((x - xo), 2)+(pow((y - yo), 2)), 0.5))
#define ST (-(x - xo)/RAD)

#define CYLINDER (1. - sqrt(sq(x - xo) + sq(y - 5.)))
double yo = 10., xo = 12.1;
int iteration = 0;

face vector muc[];
int main() {
    L0 = 25;
    mu = muc;
    NITERMAX = 100;
    /**
       For each run, the vortex structure is placed 5 units away from
       the cylinder center. The `mask`ed run.
    */
    run();
}
/**
##Implementation

for the mask method, we need to declare a boundary internal domain
(`bid`);
*/
bid cylinder;
p[cylinder] = dirichlet (1.);
//u.t[cylinder] = dirichlet (0.);
//u.n[cylinder] = dirichlet (0.);
p[left] = dirichlet(0);
p[right] = dirichlet(0);
p[bottom] = dirichlet(0);
p[top] = dirichlet(0);

pf[left] = dirichlet(0);
pf[right] = dirichlet(0);
pf[bottom] = dirichlet(0);
pf[top] = dirichlet(0);
/**
We initialize a vortex dipole with centered location `{xo, yo}` 
*/

event init (t = 0) {
    scalar psi[];
    double k = 3.83170597;
    refine (RAD < 2.0 && level <= 9);
    refine (RAD < 1.0 && level <= 10);
    refine (fabs(CYLINDER) < 0.2 && level < 9);
    refine (fabs(CYLINDER) < 0.1 && level <= 10);
    foreach()
    psi[] = ((RAD > 1)*ST/RAD +
             (RAD < 1)*(-2*j1(k*RAD)*ST/(k*j0(k)) +
                        RAD*ST));
    boundary ({psi});
    foreach() {
        u.x[] = 0.;//-((psi[0, 1] - psi[0, -1])/(2*Delta));
        u.y[] = 0.0;//(psi[1, 0] - psi[-1, 0])/(2*Delta);
    }
    /**
       For the mask method, we use the `mask()' function to mask
       cells. Note that we need to unrefine the grid (i.e. implement a
       lego boundary) to get it to run.
    */

    unrefine (y < 7.5 && level > 7);
    mask (CYLINDER > 0. ? cylinder : none);
    boundary (all);
    boundary({p, u});
    event("vtk_file");
    printf("maxiteration = %d", NITERMAX);
}
/**
We use a constant viscosity in the flow domain. This event is
compatible with all methods.
 */
event properties (i++) {
    foreach_face()
    muc.x[] = fm.x[]/500.;
    boundary ((scalar*){muc});
}
/**
Because of the spatio-temporal localization of our problem, grid adaptation is employed.
 */
event adapt (i++)
adapt_wavelet ({u.x, u.y}, (double[]){0.01, 0.01}, 10);
/**
##Output

First, movies are generated. 
*/
event movie ( t += 0.1 ; t <= 10){
scalar omega[];
vorticity (u, omega);
foreach() {
    if (x > xo)
        omega[] = level - 5;
}
output_ppm (omega, n = 512, file = "movie_cyl.mp4", min = -5, max = 5);
}

event stop (t = 10) {
}

//Output
#include "output_fields/output_vtu_foreach.h"
event vtk_file (t += 0.01; t<=10)
{
    int nf = iteration;
    scalar l[];
    foreach()
    l[] = level;

    char name[80], subname[80];
    FILE *fp;
    sprintf(name, "hs_%4.4d_n%3.3d.vtu", nf, pid());
    fp = fopen(name, "w");

    output_vtu_bin_foreach((scalar *) {l, p, pf}, (vector *) {u, muc}, 64, fp, false);
    fclose(fp);
    //  @if _MPI
    //  if (pid() == 0) {
    sprintf(name, "hs_%4.4d.pvtu", nf);
    sprintf(subname, "hs_%4.4d", nf);
    fp = fopen(name, "w");
    output_pvtu_bin((scalar *) {l, p, pf}, (vector *) {u, muc}, 64, fp, subname);
    fclose(fp);
    //  }
    //  MPI_Barrier(MPI_COMM_WORLD);
    //  @endif
    fprintf (ferr, "iteration: %d\n", iteration); fflush (ferr);
    iteration++;

}