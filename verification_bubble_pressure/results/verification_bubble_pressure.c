//#include "../src_local/centered-weugene.h"
#include "navier-stokes/centered.h"
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
// #include "navier-stokes/perfs.h"
#include "two-phase.h"
//#include "navier-stokes/conserving.h" //???should add?
#include "tension.h"

#define MAXLEVEL 10
#define MINLEVEL 4
#define feps 1e-3
#define ueps 1e-2
#define peps 1e-2
#define EPS_MAXA 1                                   // method of eps calculation

#define Re 1 // rhol*U*L/mul
#define We 1//rhol*L*U^2/sigma

#define Rrhog 10
#define Rrhos 1
#define Rmug 100
#define Rmus 1
#define dPl 0
#define dPr 0
#define Rad (0.05*L0)

#define ADAPT_SCALARS {f, p}
#define ADAPT_EPS_SCALARS {feps, peps}

p[left]   = dirichlet(0.);
pf[right] = dirichlet(0.);
p[bottom] = dirichlet(0.);
pf[top]   = dirichlet(0.);

int main() {
    L0 = 1.;
    origin (-L0/2, -L0/2.);
    N = 512;
    CFL = 0.4;
    DT = 1e-5;
    stokes = true;
    TOLERANCE = 1e-8;
    fprintf(stderr, "TOLERANCE = %g", TOLERANCE);
    rho1 = 1.0; rho2 = 1.0/Rrhog;
    mu1 = 1.0/Re;  mu2 = 1.0/(Re*Rrhog);
//  surface tension
    f.sigma = 1.0/We;
#if TREE
    f.refine = f.prolongation = fraction_refine;
#endif
    run();
}
#define U_BC 0

//p[top] = dirichlet(0);

double solution_P (double x, double y){
    return f.sigma*(sq(x) + sq(y) < sq(Rad))*(dimension - 1)/(Rad);
}

event init (t = 0) {

    if (!restore (file = "restart")) {
        int iter = 0;
        do {
            iter++;
            foreach(){
                f[] = (sq(x) + sq(y) < sq(Rad)) ? 1 : 0;
                u.x[] = 0;
//                p[] = f[]*f.sigma/Rad;
            }
            boundary ({f, u, p});
        }while (adapt_wavelet({f}, (double []){feps},
                maxlevel = MAXLEVEL, minlevel=MINLEVEL).nf != 0 && iter <= 15);
        fprintf(stderr, "init refinement iter=%d", iter);
    }else{
        fprintf(stderr, "RESTART from file");
    }
    event ("vtk_file");
}


//Output
#include "../src_local/output_vtu_foreach.h"
event end_timestep (t += 0.01){
    char subname[80]; sprintf(subname, "hrhs");
    scalar l[]; foreach() l[] = level;
    scalar Psol[]; foreach() Psol[] = solution_P(x, y);
    scalar err[]; foreach() err[] = fabs(p[] - Psol[]); //be careful with kappa, mu. They can be const unity
    output_vtu_MPI( (scalar *) {l, f, rho, p, err, pf, Psol}, (vector *) {u}, subname);
}

#if DUMP
event snapshot (i += 5000){
  char name[80];
  sprintf (name, "dump-%d", i);
  dump (file = name);
}
#endif


event adapt (i++) {
    double eps_arr[] = ADAPT_EPS_SCALARS;
//    MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) ADAPT_SCALARS,
                   eps_arr, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
}

event stop(t=100);

