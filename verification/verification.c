#include "../src_local/centered-weugene.h"
// #include "navier-stokes/perfs.h"
#define a_mu 2
#define muf1(alpha_doc, T) mu1*exp(pow(2*a_mu*y/L0,2))
#include "../src_local/rheology_model.h"
#include "tension.h"

#define MAXLEVEL 10
#define MINLEVEL 4
#define feps 1e-2
#define Teps 1e-2
#define aeps 1e-2
#define ueps 1e-3
#define mueps 1
#define EPS_MAXA 1                                   // method of eps calculation

#define Re 1 // rhol*U*L/mul
#define We 0.01//rhol*L*U^2/sigma
#define Pe 1 // Cpl*rhol*U*L/kappal
#define Ex 1 // H_tr/(Cpl*T)
#define Po 0 // A*L/U //1000
#define Aralpha 5 // E_a/RT
#define Areta 0.1 // E_eta/RT
#define CHI 26.89
#define Rrhog 10
#define Rrhos 1
#define Rmug 100
#define Rmus 1
#define Rkappag 1000
#define Rkappas 1
#define RCpg 1
#define RCps 1
#define RT 1
#define dPl 1
#define dPr 0

#define TMIN 1
#define TMAX TMIN*RT

#define ADAPT_SCALARS {mus, u.x}
#define ADAPT_EPS_SCALARS {mueps, ueps}
//#define ADAPT_SCALARS {mu.x, u.x}
//#define ADAPT_EPS_SCALARS {mueps, ueps}

int main() {
    L0 = 1.;
    origin (-L0/2, -L0/2.);
    N = 512;
    CFL = 0.4;
    DT = 1e-4;
    stokes = true;
    rho1 = 1.0; rho2 = 1.0/Rrhog; rho3 = 1.0/Rrhos;
    mu1 = 1.0/Re;  mu2 = 1.0/(Re*Rrhog);  mu3 = 1.0/(Re*Rrhos);
//  surface tension
    f.sigma = 1.0/We;
//  heat transfer model
    kappa1 = 1.0/Pe;  kappa2 = 1.0/(Pe*Rkappag);  kappa3 = 1.0/(Pe*Rkappas);
    Cp1 = 1.0;  Cp2 = 1.0/RCpg;  Cp3 = 1.0/RCps;
    Htr = Ex;
//  polymerization model
    Arrhenius_const = Po;
    Ea_by_R = Aralpha;
//  rheology model
    n_degree = 1.94;
    m_degree = 0.;
    Eeta_by_Rg = Areta;
    chi = CHI;
#if TREE
    for (scalar s in {f, fs})
    s.refine = s.prolongation = fraction_refine;
#endif
    run();
}
#define U_BC 0
#define T_BC (0.5*(TMAX + TMIN) + 0.5*(TMAX - TMIN)*tanh((x)/(L0/10.0)))

u.n[left] = neumann(0);
p[left]    = dirichlet(dPl);
pf[left]   = dirichlet(dPl);
f[left]    = dirichlet(1);
T[left]    = dirichlet(TMIN);
f[left]    = dirichlet(1);
fs[left]   = dirichlet(0);
alpha_doc[left] = dirichlet(0);

u.n[right] = neumann(0.);
p[right]   = dirichlet(dPr);
pf[right]  = dirichlet(dPr);
T[right]    = neumann(0);
f[right]    = dirichlet(1);
fs[right]   = dirichlet(0);
alpha_doc[right] = neumann(0);

u.n[top] = dirichlet(0);
u.t[top] = dirichlet(0);
alpha_doc[top] = neumann(0);
T[top] = dirichlet(T_BC);
fs[top] = neumann(0);
f[top] = neumann(0);
p[top] = neumann(0);

u.n[bottom] = dirichlet(0);
u.t[bottom] = dirichlet(0);
alpha_doc[bottom] = neumann(0);
T[bottom] = dirichlet(T_BC);
fs[bottom] = neumann(0);
f[bottom] = neumann(0);
p[bottom] = neumann(0);
double solution_u (double x, double y){
    return ((dPr-dPl)/(8*mu1*L0))*pow(L0/a_mu, 2)*(exp(-pow(a_mu,2)) - exp(-pow(-2*a_mu*y/L0,2)));
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
                T[] = TMIN;
                f[] = 1;
                fs[] = 0.0;
                sa[] = f[] + fs[];
                u.x[] = U_BC*(1.0 - fs[]);
            }
            boundary ({sa, f, fs, T, u});
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
        foreach() {
            alpha_doc[] = 0;
        }
        boundary ({alpha_doc});
    }else{
        fprintf(stderr, "RESTART from file");
    }

    event ("vtk_file");
}


//Output
#include "../src_local/output_vtu_foreach.h"
event vtk_file (i += 100){
    char subname[80]; sprintf(subname, "hrhs");
    scalar l[]; foreach() l[] = level;
    scalar usol[]; foreach() usol[] = solution_u(x, y);
    scalar err[]; foreach() err[] = fabs(u.x[] - usol[]); //be careful with kappa, mu. They can be const unity
    output_vtu_MPI( (scalar *) {l, sa, rho, p, T, alpha_doc, err, usol}, (vector *) {u, mu}, subname);
}

#if DUMP
event snapshot (i += 5000){
  char name[80];
  sprintf (name, "dump-%d", i);
  dump (file = name);
}
#endif


event adapt (i+=1000) {
    foreach() mus[] = mu.y[];
    double eps_arr[] = ADAPT_EPS_SCALARS;
//    MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) ADAPT_SCALARS,
                   eps_arr, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
}

event stop(t=200);

