#define BRINKMAN_PENALIZATION
#define DEBUG_BRINKMAN_PENALIZATION
//#define DEBUG_MINMAXVALUES
//#define DEBUG_OUTPUT_VTU_MPI

scalar fs[], omega[];
vector Us[];

#include "../src_local/centered-weugene.h"
#define RAD pow(sq(x - xo)+sq(y - yo), 0.5)
#define ST (-(x - xo)/RAD)
#define MAXLEVEL 10
#define CYLINDER (sq(x - xo) + sq(y - 5.)<1) ? 1 : 0

double yo = 10., xo = 10;
int iteration = 0;


u.n[left] = neumann(0);
p[left]    = dirichlet(37.5);
pf[left]   = dirichlet(1);
fs[left]   = dirichlet(0);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0);
pf[right]  = dirichlet(0);
fs[right]   = dirichlet(0);

u.n[top] = dirichlet(0);
u.t[top] = dirichlet(0);
fs[top] = neumann(0);
p[top] = neumann(0);

u.n[bottom] = dirichlet(0);
u.t[bottom] = dirichlet(0);
fs[bottom] = dirichlet(0);
p[bottom] = neumann(0);
face vector muv[];
int main() {
    L0 = 1;
    origin(-L0/2,-L0/2);
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
            fs[] = (sq(x) +sq(y) < sq(0.1*L0));
            foreach_dimension() {Us.x[] = 0;}
            u.x[] = 1*(1-fs[]); u.y[] = 0;
            p[] = 0;
        }
        boundary ({fs, Us, u.x, u.y, p});
        if (it>=10) printf("WARNING: does not converge... ");
    }while (adapt_wavelet({fs, p}, (double []){0.01, 0.01, 0.01},
                          maxlevel = MAXLEVEL, minlevel = 3).nf != 0 && it <= 10);
    event ("vtk_file");
}

event properties (i++)
{
    foreach_face()
    muv.x[] = 1;
}

#define ADAPT_SCALARS {fs, omega, p}
#define ADAPT_EPS_SCALARS {0.01, 0.01, 0.01}
//Because of the spatio-temporal localization of our problem, grid adaptation is employed.
event adapt (i++){
    double eps_arr[] = ADAPT_EPS_SCALARS;
    MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = MAXLEVEL, minlevel = 4);
}
//Output
#include "../src_local/output_vtu_foreach.h"
event end_timestep (i++){
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
