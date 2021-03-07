#include "../src_local/centered-weugene.h"
//#include "navier-stokes/centered.h"
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
// #include "navier-stokes/perfs.h"
#include "two-phase.h"
//#include "navier-stokes/conserving.h" //???should add?
#define NOT_CFL_SIGMA
#include "../src_local/tension.h"

#define MINLEVEL 4
#define EPS_MAXA 1                                   // method of eps calculation
#define Rad (0.25)
#define tend 3

#define ADAPT_SCALARS {gf.x, gf.y, fwall}
#define ADAPT_EPS_SCALARS {feps, feps, feps}

int MAXLEVEL = 9;
double feps = 1e-2, ueps = 1e-2, peps = 1e-2;
double grav = 0.98;
vector gf[];
scalar fwall[];
scalar omega[];

p[bottom] = dirichlet(0);
u.n[bottom] = dirichlet(0);
u.t[bottom] = dirichlet(0);
u.n[top] = dirichlet(0.);
u.t[top] = dirichlet(0.);

int main(int argc, char * argv[]) {
    double eps=1e-2;
    if (argc > 1) {
        eps = atoi(argv[1]); //convert from string to int
    }
    feps = eps; ueps = eps; peps = eps;
    L0 = 2.;
    origin (-0.5, 0);
    N = 512;
    CFL = 0.4;
    DT = 3e-4;
    stokes = true;
    TOLERANCE = 1e-8;
    fprintf(stderr, "TOLERANCE = %g", TOLERANCE);
#if TREE
    f.refine = f.prolongation = fraction_refine;
    fwall.refine = fwall.prolongation = fraction_refine;
#endif
////CASE A
//    rho1 = 1000; rho2 = 100;
//    mu1 = 10;  mu2 = 1;
//    f.sigma = 24.5;
//    run();
//CASE B
    rho1 = 1000; rho2 = 1;
    mu1 = 10;  mu2 = 0.1;
    f.sigma = 1.96;
    run();
}

event init (t = 0) {

    if (!restore (file = "restart")) {
        int iter = 0;
        do {
            iter++;
            foreach(){
                f[] = (sq(x-0.5) + sq(y-0.5) > sq(Rad)) ? 1 : 0;
                fwall[] = (x<0 || x>1) ? 1 : 0;
                u.x[] = 0;
            }
            boundary ({f, u});
            gradients ({f}, {gf});
        }while (adapt_wavelet({gf.x, gf.y, fwall}, (double []){feps, feps, feps},
                maxlevel = MAXLEVEL, minlevel=MINLEVEL).nf != 0 && iter <= 15);
        fprintf(stderr, "init refinement iter=%d", iter);
    }else{
        fprintf(stderr, "RESTART from file");
    }
    event ("vtk_file");
}

// Gravity
event acceleration(i++){
face vector av = a;
foreach_face(y){
        av.y[] -= grav;
}
}
event velocity_correction(i++){
    foreach() foreach_dimension() u.x[] *= (1.0 - fwall[]);
    boundary({u});
}
//Output
#include "../src_local/output_vtu_foreach.h"
event end_timestep (t += 0.01){
    char subname[80]; sprintf(subname, "hrhs");
    scalar l[]; foreach() l[] = level;
    vorticity (u, omega);
    output_vtu_MPI( (scalar *) {l, f, fwall, omega, rho, p}, (vector *) {u, a}, subname);
}

#if DUMP
event snapshot (i += 5000){
  char name[80];
  sprintf (name, "dump-%d", i);
  dump (file = name);
}
#endif


event adapt (i++) {
    gradients ({f}, {gf});
    double eps_arr[] = ADAPT_EPS_SCALARS;
    MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) ADAPT_SCALARS,
                   eps_arr, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
}

event stop(t=tend);

//                f[] = (sq(x) + sq(y) > sq(Rad)) && (y<0.2*L0) ? 1 : 0;