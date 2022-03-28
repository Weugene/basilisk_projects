#include "../src_local/centered-weugene.h"
// #include "navier-stokes/perfs.h"
#include "../src_local/not-used/three-phase-weugene.h"
#include "tension.h"

#define MAXLEVEL 10
#define MINLEVEL 4
#define feps 1e-2
#define ueps 1e-2
#define EPS_MAXA 1                                   // method of eps calculation

#define Re 1 // rhol*U*L/mul
#define We 1//rhol*L*U^2/sigma

#define Rrhog 10
#define Rrhos 1
#define Rmug 100
#define Rmus 1
#define dPl 0
#define dPr 0


#define ADAPT_SCALARS {f, u.x}
#define ADAPT_EPS_SCALARS {feps, ueps}

int main() {
    L0 = 1.;
    origin (-L0/2, -L0/2.);
    N = 512;
    CFL = 0.4;
    DT = 1e-5;
    stokes = true;
    rho1 = 1.0; rho2 = 1.0/Rrhog; rho3 = 1.0/Rrhos;
    mu1 = 1.0/Re;  mu2 = 1.0/(Re*Rrhog);  mu3 = 1.0/(Re*Rrhos);
//  surface tension
    f.sigma = 1.0/We;
#if TREE
    for (scalar s in {f, fs})
    s.refine = s.prolongation = fraction_refine;
#endif
    run();
}
#define U_BC 0

p[top] = dirichlet(0);

double solution_P (double x, double y){
    return f.sigma*(sq(x) + sq(y) < sq(L0/4))*(dimension - 1)/(L0/4);
}

scalar sa[];
vector gva[];
scalar mus[];

event init (t = 0) {
    double eps_arr[] = ADAPT_EPS_SCALARS;
    size_t n_eps_arr = sizeof(eps_arr) / sizeof(double);
    if (!restore (file = "restart")) {
        int iter = 0;
        do {
            iter++;
            foreach(){
                f[] = (sq(x) + sq(y) < sq(L0/4)) ? 1 : 0;
                fs[] = 0.0;
                sa[] = f[] + fs[];
                u.x[] = U_BC*(1.0 - fs[]);
            }
            boundary ({sa, f, fs, u});
            gradients ({sa}, {gva});
            event("properties");
            foreach() mus[] = mu.y[];
            double eps_arr_tmp[] = ADAPT_EPS_SCALARS;
            MinMaxValues(ADAPT_SCALARS, eps_arr_tmp);
            for (int k=0; k<n_eps_arr; k++)
                eps_arr[k] = eps_arr_tmp[k];
        }while (adapt_wavelet((scalar *)ADAPT_SCALARS, eps_arr,
                maxlevel = MAXLEVEL, minlevel=MINLEVEL).nf != 0 && iter <= 15);
        fprintf(stderr, "init refinement iter=%d", iter);
    }else{
        fprintf(stderr, "RESTART from file");
    }

    event ("vtk_file");
}


//Output
#include "../src_local/output_vtu_foreach.h"
event vtk_file (i += 1){
    char subname[80]; sprintf(subname, "hrhs");
    scalar l[]; foreach() l[] = level;
    scalar Psol[]; foreach() Psol[] = solution_P(x, y);
    scalar err[]; foreach() err[] = fabs(p[] - Psol[]); //be careful with kappa, mu. They can be const unity
    output_vtu_MPI( (scalar *) {l, rho, p, err, pf, Psol}, (vector *) {u}, subname);
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

