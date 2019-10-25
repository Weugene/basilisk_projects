//#include "../src_local/centered-weugene.h"
#include "navier-stokes/centered.h"
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
// #include "navier-stokes/perfs.h"
#include "two-phase.h"
//#include "navier-stokes/conserving.h" //???should add?
#include "tension.h"
#define MINLEVEL 4

#define EPS_MAXA 1                                   // method of eps calculation
#define Rad (0.1*L0)
#define r_d (0.1*Rad)

#define ADAPT_SCALARS {f, p}
#define ADAPT_EPS_SCALARS {feps, peps}

int MAXLEVEL = 10;
double feps = 1e-3, ueps = 1e-2, peps = 1e-2;
double Re = 1, We = 1, Bo=1; // Re = rhol*U*L/mul //We = rhol*L*U^2/sigma
double Rrhog = 10, Rmug = 100;
p[left]   = dirichlet(0.);
pf[right] = dirichlet(0.);
p[bottom] = dirichlet(0.);
pf[top]   = dirichlet(0.);

int main(int argc, char * argv[]) {
//    maxruntime (&argc, argv);
    if (argc > 1) {
        MAXLEVEL = atoi(argv[1]); //convert from string to int
    }else if (argc > 2) {
        feps = atof(argv[2]);
    }else if (argc > 3) {
        ueps = atof(argv[3]);
    }else if (argc > 4) {
        peps = atof(argv[4]);
    }else if (argc > 5) {
        Re = atof(argv[5]);
    }else if (argc > 6) {
        We = atof(argv[6]);
    }else if (argc > 7) {
        Rrhog = atof(argv[7]);
    }else if (argc > 8) {
        Rmug = atof(argv[8]);
    }
    L0 = 1.;
    origin (-L0/2, -L0/2.);
    N = 512;
    CFL = 0.4;
    DT = 1e-5;
    stokes = true;
    TOLERANCE = 1e-8;
    fprintf(stderr, "TOLERANCE = %g", TOLERANCE);
    rho1 = 1.0; rho2 = 1.0/Rrhog;
    mu1 = 1.0/Re;  mu2 = 1.0/(Re*Rmug);
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
        double r2, r2i, theta;
        do {
            iter++;
            foreach(){
                theta = atan(y/x);
                r2 = (sq(x) + sq(y));
                r2i = sq(Rad) * (1 + r_d * cos(10 * theta));
                f[] = (r2 > r2i)? 1 : 0;
                u.x[] = 0;
                pf[] = f[] * f.sigma/Rad;//??
            }
            boundary ({f, u, pf});
        }while (adapt_wavelet({f}, (double []){feps},
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
        av.y[] -= Bo;
}
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

