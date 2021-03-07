//A test for curved-boundary implementations
//On this page a comparison between curved-boundary-implementations is presented.
// Considering mask, Stephane’s trick and the embedded boundary method.
// The test case considers a two-dimensional vortex-cylinder collision,
// and we consider the evolution of the enstrophy as a critical statistic for the representation of the flow.

#include "embed.h" //The volume and area fractions are stored in these fields:cs, fs
#include "navier-stokes/centered.h"
#include "navier-stokes/double-projection.h"
#define RAD pow(sq(x - xo)+sq(y - yo), 0.5)
#define ST (-(x - xo)/RAD)
#define MAXLEVEL 10
#define CYLINDER (1. - sqrt(sq(x - xo) + sq(y - 5.)))* \
                 (1. - sqrt(sq(x - xo-3) + sq(y - 5.)))* \
                 (1. - sqrt(sq(x - xo+3) + sq(y - 5.)))* \
                 (1. - sqrt(sq(x - xo-6) + sq(y - 5.)))* \
                 (1. - sqrt(sq(x - xo+6) + sq(y - 5.)))
double yo = 20., xo = 12.1;
int j;

face vector muc[];
int main() {
    L0 = 25;
    mu = muc;
//    For each run, the vortex structure is placed 5 units away from the cylinder center.
//    The masked run has label j = 0.
    j = 0;
//    run();
//    Next, we try “Stephane’s trick”.
    j++;
    run();
//    Finally, there is the embedded method.
    j++;
    run();
}
//Implementation
//for the mask method, we need to declare a boundary internal domain (bid);

//bid cylinder;
//u.t[cylinder] = dirichlet (0.);
//We initialize a vortex dipole with centered location {xo, yo}

event init (t = 0) {
    scalar psi[];
    double k = 3.83170597;
    refine (RAD < 2.0 && level <= 9);
    refine (RAD < 1.0 && level <= 10);
    refine (fabs(CYLINDER) < 0.2 && level < 9); //in the vicinity cylinder
    refine (fabs(CYLINDER) < 0.1 && level <= 10);
    foreach()
    psi[] = ((RAD > 1)*ST/RAD +
             (RAD < 1)*(-2*j1(k*RAD)*ST/(k*j0(k)) + //j0, j1 Bessel functions
                        RAD*ST));
    boundary ({psi});
    foreach() {
        u.x[] = -((psi[0, 1] - psi[0, -1])/(2*Delta)); //ux = -dyPsi
        u.y[] = (psi[1, 0] - psi[-1, 0])/(2*Delta);    //uy =  dxPsi
    }
//    For the mask method, we use the `mask()’ function to mask cells.
//    Note that we need to unrefine the grid (i.e. implement a lego boundary) to get it to run.

//    if (j == 0) {
//        unrefine (y < 7.5 && level > 7);
//        mask (CYLINDER > 0. ? cylinder : none);
//    }
//    For the embedded boundary, we compute its location and implement it like so:

    if (j == 2) {
        scalar phi[];
        foreach_vertex()
        phi[] = -CYLINDER; // phi<0 in, phi>0 out
        fractions (phi, cs, fs);//=>cs and fs are calculated
        u.n[embed] = dirichlet (0.);
        u.t[embed] = dirichlet (0.);
    }
    boundary (all);
}
//We use a constant viscosity in the flow domain. This event is compatible with all methods.

event properties (i++) {
    foreach_face()
    {
        muc.x[] = fm.x[] / 500.;//0.002=1/500
    }
    boundary ((scalar*){muc});
}
//Stephane’s trick is implemented via an additional event.

event Stephanes_trick (i++) {
    if (j == 1) {
        scalar f[];
        fraction (f, CYLINDER);//f=0 out. f=1 in

        foreach(){
            foreach_dimension()
            u.x[] -= u.x[]*f[];//damping of velocity
        }
    }
}
//Because of the spatio-temporal localization of our problem, grid adaptation is employed.

event adapt (i++){
    adapt_wavelet ({cs, u.x, u.y}, (double[]){0.001, 0.001, 0.001}, MAXLEVEL);
}
//Output
//        First, movies are generated.

event movie ( t += 0.1 ; t <= 10){
scalar omega[], l;
vorticity (u, omega);
foreach() {
    omega[] = (x > xo)?level:fabs(omega[]);
}
output_ppm (omega, n = 1024, file = "movie_cyl.mp4", min = 0, max = MAXLEVEL); // n=512 number of pixels
}
//The overall dynamics appear quite similar.
//The movie cycles over the three methods
//Second, we quantify the total vorticity via the enstrophy (E),

event diag (i += 5) {
    double E = 0;
    boundary ({u.x, u.y});
    scalar omega[];
    vorticity (u , omega);
    foreach(){
        double vort = omega[];
        if (cs[] < 1. && cs[] > 0){ //Embedded boundary cell //f=0 out. f=1 in
            coord b, n;
            double area = embed_geometry (point, &b, &n);
            vort = embed_vorticity (point, u, b, n);
        }
        E += dv()*sq(vort);
    }
    char fname[99];
    sprintf (fname, "data_cyl%d", j);
    static FILE * fp = fopen (fname, "w");
    fprintf (fp, "%d\t%g\t%g\n", i, t, E);
    fflush (fp);
}
//The maximum value is sensitive to the boundary-implementation details, the value at t = 10 appears more robust.
//The maximum value is sensitive to the boundary-implementation details, the value at t = 10 appears more robust.

//Finally we compare how long it takes for each run to complete.

event stop (t = 10) {
    static FILE * fp = fopen ("perf_cyl", "w");
    timing s = timer_timing (perf.gt, iter, perf.tnc, NULL);
    fprintf (fp, "%d\t%g\t%d\t%g\n", j, s.real, i, s.speed);
    fflush (fp);
    return 1;
}
