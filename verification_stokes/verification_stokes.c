#define BRINKMAN_PENALIZATION 1
#define DEBUG_BRINKMAN_PENALIZATION
//#define DEBUG_MINMAXVALUES
//#define DEBUG_OUTPUT_VTU_MPI

scalar fs[];
vector n_solv[];
#include "../src_local/centered-weugene.h"

#define MAXLEVEL 10
#define Re 1
#define U_inf (Re/diam)
#define Delta_int (rad/100.)
#undef SEPS
#define SEPS 1e-30
#define JACOBI 1
double yo = 0, xo = 0, rad = 0.5, diam = 1;
int iteration = 0;
scalar omega[];
face vector muv[];

u.n[left] = dirichlet(U_inf);
u.t[left] = dirichlet(0);

u.n[right] = neumann(0);
u.t[right] = neumann(0);

u.n[bottom] = neumann(0);
u.t[bottom] = neumann(0);

u.n[top] = neumann(0);
u.t[top] = neumann(0);

p[left]   = neumann(0);
p[right]  = dirichlet(0);
p[top]    = neumann(0);
p[bottom] = neumann(0);

int main() {
    L0 = 10;
    origin(-L0/3, -L0/2);
    eta_s = 1e-3;
    nu_s = 0;
    TOLERANCE = 1e-5;
    DT = 1e-4;
    N = 512;
    mu = fm;
//    mu = muv;
    n_sol = n_solv;
    stokes = true;
    run();
}

event properties (i++) {
//    foreach_face() {
//        muv.x[] = fm.x[];
//    }
    foreach() {
        theta = atan2(y, x + SEPS);
        n_solv.x[] = -cos(theta);
        n_solv.y[] = -sin(theta);
    }
}

event init (t = 0) {
    int it = 0;
    double theta, r2;
    do {
        it++;
        foreach() {
            r2 = sq(x - xo) + sq(y - yo);
            theta = atan2(y, x + SEPS);
            fs[] = 0.5*(1 - tanh((r2- sq(rad))/Delta_int));
            u.x[] = U_inf*(1.0 - fs[]);
            u.y[] = 0;
        }
        boundary ({fs, u.x, u.y});
        if (it>=10) printf("WARNING: does not converge... ");
    }while (adapt_wavelet({fs, u.x, u.y}, (double []){1e-3, 1e-3, 1e-3},
                          maxlevel = MAXLEVEL, minlevel = 3).nf != 0 && it <= 10);

    event ("properties");
    event ("end_timestep");
}

void exact_solution(vector ve, scalar pe){
    double theta, r, vr, vth;
    foreach() {
        r = sqrt(sq(x) + sq(y));
        theta = atan2(y, x + SEPS);
      #if 0
        vr  =   (r > rad) * U_inf * (1 - sq(rad/(r + SEPS))) * cos(theta);
        vth = - (r > rad) * U_inf * (1 + sq(rad/(r + SEPS))) * sin(theta);
        ve.x[] = vr * cos(theta) - vth * sin(theta);
        ve.y[] = vr * sin(theta) + vth * cos(theta);
        pe[] = 0.5 * rho[] * (sq(U_inf) - sq(ve.x[]) - sq(ve.y[]));
      #else
        ve.x[] = (r > rad) ? ((sq(R) - sq(r))*sq(cos(theta)) + sq(r)*log(r/R) + 0.5*(sq(r) - sq(R)))/sq(r) : 0;
        ve.y[] = (r > rad) ? (sq(R) - sq(r))*cos(theta)*sin(theta)/sq(r) : 0;
        pe[]   = - 2*cos(theta)/r;
      #endif

    }
}

//Output
#include "../src_local/output_vtu_foreach.h"
event end_timestep (i += 10){
    char subname[80]; sprintf(subname, "br");
    scalar l[], pe[];
    vector ve[];
    exact_solution(ve, pe);
    vorticity (u, omega);
    foreach() l[] = level;
    output_vtu_MPI( (scalar *) {l, omega, fs, p, pe}, (vector *) {u, total_rhs, dbp, ve, n_sol}, subname, 0);
}

#define ADAPT_SCALARS {fs, omega}
#define ADAPT_EPS_SCALARS {1e-3, 1e-2}
//Because of the spatio-temporal localization of our problem, grid adaptation is employed.
event adapt (i++){
    double eps_arr[] = ADAPT_EPS_SCALARS;
    MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = MAXLEVEL, minlevel = 4);
    boundary ({u.x, u.y});
}

event log_file(i += 10){
    fprintf(ferr, "i=%d t=%g dt=%g ui=%d un=%d pi=%d pn=%d pfi=%d pfn=%d ",
            i, t, dt, mgu.i, mgu.nrelax, mgp.i, mgp.nrelax, mgpf.i, mgpf.nrelax);
}
event stop(i = 1000);
