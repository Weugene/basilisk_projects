#include "navier-stokes/centered.h"
#include "navier-stokes/double-projection.h"
#include "fractions.h"
#define RAD pow(sq(x - xo)+sq(y - yo), 0.5)
#define ST (-(x - xo)/RAD)
#define MAXLEVEL 10
#define CYLINDER (1. - sqrt(sq(x - xo) + sq(y - 5.)))
//* \
//                 (1. - sq(sq(x - xo-3) + sq(y - 5.)))* \
//                 (1. - sq(sq(x - xo+3) + sq(y - 5.)))* \
//                 (1. - sq(sq(x - xo-6) + sq(y - 5.)))* \
//                 (1. - sq(sq(x - xo+6) + sq(y - 5.)))
double yo = 10., xo = 12.5;
int iteration = 0;
scalar channel[];
face vector channelf[];
int main() {
    L0 = 25;
    run();
}

event init (t = 0) {
    scalar psi[];
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
        // fraction of Obstacle 1 in cylinder, 0 outside


        vertex scalar phi[];
        foreach_vertex()
          phi[] = CYLINDER;
        boundary ({phi});
        fractions (phi, channel, channelf);
        if (it>=10) printf("WARNING: does not converge... ");
    }while (adapt_wavelet({channel, u.x, u.y}, (double []){0.001, 0.01, 0.01},
                          maxlevel = MAXLEVEL, minlevel = 3).nf != 0 && it <= 10);

    boundary (all);
    event ("vtk_file");
}

//event bc (i += 1) {
//    //no slip boundary conditions
//    foreach()
//        foreach_dimension()
//            u.x[] = (1 - channel[])*u.x[];
//    boundary ((scalar *){u , p});
//}

const double eta = 0.01;
//event acceleration (i++) {
//    face vector av = a; //pointer equals?
//
//    foreach_face(x)
//    av.x[] = 1;//-channelf.x[] * uf.x[] / eta;
//
////    foreach_face(y)
////    av.y[] = -channel[] * u.y[] / eta;
////#if dimension==3
////    foreach_face(z)
////        av.z[] = -channel[] * u.z[] / eta;
////#endif
//
//    boundary ((scalar *){av});
//}

//Because of the spatio-temporal localization of our problem, grid adaptation is employed.
event adapt (i++){
    adapt_wavelet ({channel, u.x, u.y}, (double[]){0.001, 0.01, 0.01}, MAXLEVEL);
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

    output_vtu_bin_foreach((scalar *) {l, omega, channel}, (vector *) {u}, 64, fp, false);
    fclose(fp);
    @if _MPI
    if (pid() == 0) {
        sprintf(name, "hs_%4.4d.pvtu", nf);
        sprintf(subname, "hs_%4.4d", nf);
        fp = fopen(name, "w");
        output_pvtu_bin((scalar *) {l, omega, channel}, (vector *) {u}, 64, fp, subname);
        fclose(fp);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    @endif

    iteration++;
}
