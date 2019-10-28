#define BRINKMAN_PENALIZATION
#define DEBUG_BRINKMAN_PENALIZATION
//#define DEBUG_MINMAXVALUES
//#define DEBUG_OUTPUT_VTU_MPI

scalar fs[], omega[];
vector Us[];

#include "../src_local/centered-weugene.h"

#define MAXLEVEL 10
#define Re 200
#define CYLINDER (sq(x - xo) + sq(y - yo)< sq(rad)) ? 1 : 0

double yo = 1.65, xo = 1.5, rad = 0.05;
int iteration = 0;


u.n[bottom] = dirichlet(1);
u.n[top] = neumann(0.);
p[top]   = dirichlet(0);


face vector muv[];
int main() {
    L0 = 3;
//    origin(-L0/2,-L0/2);
    eta_s =1e-5;
    TOLERANCE = 1e-5;
    DT = 1e-4;
    mu = muv;
    run();
}

event init (t = 0) {
    int it = 0;
    do {
        it++;
        // fraction of Obstacle 1 in cylinder, 0 outside
        foreach() {
            fs[] = CYLINDER;
            foreach_dimension() {Us.x[] = 0;}
            u.y[] = (1 + 0.01*noise())*(1-fs[]); u.x[] = 0;
        }
        boundary ({fs, Us, u.x, u.y});
        if (it>=10) printf("WARNING: does not converge... ");
    }while (adapt_wavelet({fs, u.x, u.y}, (double []){0.01, 0.01, 0.01},
                          maxlevel = MAXLEVEL, minlevel = 3).nf != 0 && it <= 10);
    event ("end_timestep");
}

event properties (i++)
{
    foreach_face()
    muv.x[] = 2*rad/Re;
}

#define ADAPT_SCALARS {fs, omega}
#define ADAPT_EPS_SCALARS {0.01, 0.01}
//Because of the spatio-temporal localization of our problem, grid adaptation is employed.
event adapt (i++){
    double eps_arr[] = ADAPT_EPS_SCALARS;
    MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = MAXLEVEL, minlevel = 4);
}
//Output
#include "../src_local/output_vtu_foreach.h"
event end_timestep (t += 0.1){
//event end_timestep (t += 0.01){
    char subname[80]; sprintf(subname, "br");
    scalar l[];
    vorticity (u, omega);
    foreach() l[] = level;
    output_vtu_MPI( (scalar *) {l, omega, fs, p}, (vector *) {u, dbp}, subname);
}

event stop(t = 10);


//    scalar psi[];
//    double k = 3.83170597;
//setting scalar functions for defining velocity field. ux = -dyPsi uy =  dxPsi
//        foreach()
//        psi[] = ((RAD > 1)*ST/RAD +
//                 (RAD < 1)*(-2*j1(k*RAD)*ST/(k*j0(k)) + //j0, j1 Bessel functions
//                            RAD*ST));
//        boundary ({psi});
//        foreach() {
//            u.x[] = -((psi[0, 1] - psi[0, -1])/(2*Delta)); //ux = -dyPsi
//            u.y[] = (psi[1, 0] - psi[-1, 0])/(2*Delta);    //uy =  dxPsi
//        }
//u.n[left] = neumann(0);
//p[left]    = dirichlet(37.5);
//pf[left]   = dirichlet(1);
//fs[left]   = dirichlet(0);
//
//u.n[right] = neumann(0.);
//p[right]   = dirichlet(0);
//pf[right]  = dirichlet(0);
//fs[right]   = dirichlet(0);
//
//u.n[top] = dirichlet(0);
//u.t[top] = dirichlet(0);
//fs[top] = neumann(0);
//p[top] = neumann(0);
//
//u.n[bottom] = dirichlet(0);
//u.t[bottom] = dirichlet(0);
//fs[bottom] = dirichlet(0);
//p[bottom] = neumann(0);





//vertex scalar phi[];
//int it = 0;
//do {
//it++;
//foreach_vertex() {
//    phi[] = intersection (0.5 - y, 0.5 + y);
//    phi[] = intersection (phi[], sq(x) + sq(y) - sq(0.125/2.));
//}
//boundary ({phi});
//fractions (phi, cs, fs);
//foreach() {
//    u.x[] = cs[] ? 1. : 0.;
//    u.y[] = 0;
//}
//if (it>=10) printf("WARNING: does not converge... ");
//}while (adapt_wavelet({cs, u.x, u.y}, (double []){1e-3, 1e-3, 1e-3},
//maxlevel = MAXLEVEL, minlevel = 4).nf != 0 && it <= 10);
//event ("end_timestep");