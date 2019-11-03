#define BRINKMAN_PENALIZATION
#define DEBUG_BRINKMAN_PENALIZATION
//#define DEBUG_MINMAXVALUES
//#define DEBUG_OUTPUT_VTU_MPI

scalar fs[], omega[];
vector Us[];

#include "../src_local/centered-weugene.h"

#define MAXLEVEL 10
#define Re 40
//#define U_inf Re/(diam)
#define CYLINDER (sq(x - xo) + sq(y - yo)< sq(rad)) ? 1 : 0

double yo = 16.5, xo = 15, rad = 0.5, diam = 1;
int iteration = 0;


//u.n[bottom] = dirichlet(U_inf*min(t, 1));
//u.n[top] = neumann(0.);
p[top]   = dirichlet(0);


int main() {
    L0 = 30;
    eta_s =1e-7;
    TOLERANCE = 1e-5;
    DT = 1e-4;
    const face vector g[] = {0,12.7388};
    a = g;
    stokes = true;
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
            u.y[] = 0; u.x[] = 0;
        }
        boundary ({fs, Us, u.x, u.y});
        if (it>=10) printf("WARNING: does not converge... ");
    }while (adapt_wavelet({fs, u.x, u.y}, (double []){0.001, 0.001, 0.001},
                          maxlevel = MAXLEVEL, minlevel = 3).nf != 0 && it <= 10);
    event ("end_timestep");
}



#define ADAPT_SCALARS {fs, omega}
#define ADAPT_EPS_SCALARS {0.001, 0.001}
//Because of the spatio-temporal localization of our problem, grid adaptation is employed.
event adapt (i++){
    double eps_arr[] = ADAPT_EPS_SCALARS;
    MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = MAXLEVEL, minlevel = 4);
    boundary ({u.x, u.y});
}
//Output
#include "../src_local/output_vtu_foreach.h"
event end_timestep (t += 0.01){
//event end_timestep (i += 1){
    char subname[80]; sprintf(subname, "br");
    scalar l[];
    vorticity (u, omega);
    foreach() l[] = level;
    output_vtu_MPI( (scalar *) {l, omega, fs, p}, (vector *) {u, dbp}, subname, 0);
}

event stop(t = 10);
