#include "embed.h" //The volume and area fractions are stored in these fields:cs, fs
#include "navier-stokes/centered.h"
#include "navier-stokes/double-projection.h"
#define RAD pow(sq(x - xo)+sq(y - yo), 0.5)
#define ST (-(x - xo)/RAD)
#define MAXLEVEL 10
#define CYLINDER (1. - sqrt(sq(x - xo) + sq(y - 5.)))
//* \
//                 (1. - sqrt(sq(x - xo-3) + sq(y - 5.)))* \
//                 (1. - sqrt(sq(x - xo+3) + sq(y - 5.)))* \
//                 (1. - sqrt(sq(x - xo-6) + sq(y - 5.)))* \
//                 (1. - sqrt(sq(x - xo+6) + sq(y - 5.)))
double yo = 10., xo = 12.5;
int iteration = 0;
face vector muc[];
int main() {
    L0 = 25;
    mu = muc;
    run();
}

event init (t = 0) {
    scalar psi[], phi[];
    double k = 3.83170597;

    int it = 0;

    do {
        it++;
        printf("u\n");
        //setting scalar functions for defining velocity field. ux = -dyPsi uy =  dxPsi
        foreach()
        psi[] = ((RAD > 1)*ST/RAD +
                 (RAD < 1)*(-2*j1(k*RAD)*ST/(k*j0(k)) + //j0, j1 Bessel functions
                            RAD*ST));
        boundary ({psi});
        foreach() {
            u.x[] = -((psi[0, 1] - psi[0, -1])/(2*Delta)); //ux = -dyPsi
            u.y[] = (psi[1, 0] - psi[-1, 0])/(2*Delta);    //uy =  dxPsi
        }

        printf("channel\n");
        foreach_vertex()
            phi[] = -CYLINDER; // phi<0 in, phi>0 out of cyl
        fractions (phi, cs, fs);//=>cs and fs are calculated to define an obstacle//cs, fs each direction 0 in cylinder, 1 out of cyl
        //fraction(channel, RAD - 1);
        if (it>=10) printf("WARNING: does not converge... ");
    }while (adapt_wavelet({cs, u.x, u.y}, (double []){0.001, 0.01, 0.01},
                          maxlevel = MAXLEVEL, minlevel = 3).nf != 0 && it <= 10);


//    For the embedded boundary, we compute its location and implement it like so:


    u.n[embed] = dirichlet (0.);
    u.t[embed] = dirichlet (0.);

    boundary (all);
    event ("vtk_file");
}
//We use a constant viscosity in the flow domain. This event is compatible with all methods.

event properties (i++) {
    foreach_face()
    {
        muc.x[] = fm.x[] / 500.;//0.002=1/500
    }
    boundary ((scalar*){muc});
}

//Because of the spatio-temporal localization of our problem, grid adaptation is employed.
event adapt (i++){
    adapt_wavelet ({cs, u.x, u.y}, (double[]){0.001, 0.01, 0.01}, MAXLEVEL);
}
//Output
#include "output_fields/output_vtu_foreach.h"
event vtk_file (t += 0.1; t<=100)
{
    int nf = iteration;
    scalar omega[], l[];
    vorticity (u, omega);
    foreach()
        l[] = level;

    char name[80], subname[80];
    FILE *fp;
    sprintf(name, "hs_%4.4d_n%3.3d.vtu", nf, pid());
    fp = fopen(name, "w");

    output_vtu_bin_foreach((scalar *) {l, omega, cs}, (vector *) {u, muc, fs}, 64, fp, false);
    fclose(fp);
    @if _MPI
    if (pid() == 0) {
        sprintf(name, "hs_%4.4d.pvtu", nf);
        sprintf(subname, "hs_%4.4d", nf);
        fp = fopen(name, "w");
        output_pvtu_bin((scalar *) {l, omega, cs}, (vector *) {u, muc, fs}, 64, fp, subname);
        fclose(fp);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    @endif

    iteration++;
}
