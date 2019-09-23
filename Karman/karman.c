#include "../src_local/centered-weugene.h"
// #include "navier-stokes/perfs.h"
#include "../src_local/rheology_model.h"
#include "tension.h"

#define MAXLEVEL 10
#define MINLEVEL 4
#define feps 1e-2
#define Teps 1e-2
#define aeps 1e-2
#define ueps 1e-2


#define Re 1 // rhol*U*L/mul
#define We 0.01//rhol*L*U^2/sigma
#define Pe 1 // Cpl*rhol*U*L/kappal
#define Ec 0.01 // U^2/(Cpl*T)
#define Ex 1 // H_tr/(Cpl*T)
#define Po 1000 // A*L/U //1000
#define Aralpha 5 // E_a/RT
#define Areta 0.1 // E_eta/RT
#define Rrhog 10
#define Rrhos 1
#define Rmug 100
#define Rmus 1
#define Rkappag 1000
#define Rkappas 1
#define RCpg 1
#define RCps 1
#define RT 2



#define TMIN 1
#define TMAX TMIN*RT


int main() {
    L0 = 1.;
    origin (-L0/2, -L0/2.);
    N = 512;
    CFL = 0.5;
    DT = 1e-6;
    stokes = true;
//    rho1 = 1140; rho2 = 1; rho3 = 2000;
//    mu1 = 0.155; mu2 = 1.81e-5; mu3 = 1;
//double Arrhenius_const = 4.53e+7;//1/s
//double Ea_by_R = 72900/8.314;// Kelvin
    rho1 = 1.0; rho2 = 1.0/Rrhog; rho3 = 1.0/Rrhos;
    mu1 = 1.0/Re;  mu2 = 1.0/(Re*Rrhog);  mu3 = 1.0/(Re*Rrhos);
//  surface tension
//    f.sigma = 1.0/We;
//  heat transfer model
    kappa1 = 1.0/Pe;  kappa2 = 1.0/(Pe*Rkappag);  kappa3 = 1.0/(Pe*Rkappas);
    Cp1 = 1.0;  Cp2 = 1.0/RCpg;  Cp3 = 1.0/RCps;
    Htr = Ex;
//  polymerization model
    Arrhenius_const = Po;
    Ea_by_R = Aralpha;
//  rheology model
    n_degree = 1.667;
    m_degree = 0.333;
    Eeta_by_Rg = Areta;
    chi = 10;
    for (scalar s in {f,fs}) {
        s.refine = s.prolongation = fraction_refine;
    }
    run();
}
#define U_BC (sqrt(1 - pow(y/(L0/2.),8)))
//#define UT_BC exp(-pow(x + L0/2.,2)/(2*pow(L0/10.,2)))
#define T_BC (0.5*(TMAX + TMIN) + 0.5*(TMAX - TMIN)*tanh((x)/(L0/10.0)))
u.n[left]  = dirichlet(U_BC);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(1);
//f[left]    = dirichlet((fabs(y)<sin(0.5*pi*t))?0:1);
T[left]    = dirichlet(TMIN);
fs[left]   = dirichlet(0);
alpha_doc[left] = 0;//dirichlet(0);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);
T[right]    = neumann(0);
fs[right]   = neumann(0);
alpha_doc[right] = neumann(0);


u.n[top] = dirichlet(0);
u.t[top] = dirichlet(0);
//u.t[top] = dirichlet(0);
alpha_doc[top] = neumann(0);
T[top] = dirichlet(T_BC);
fs[top] = neumann(0);
f[top] = neumann(0);

u.n[bottom] = dirichlet(0.);
u.t[bottom] = dirichlet(0);
//u.t[bottom] = dirichlet(0);
alpha_doc[bottom] = neumann(0);
T[bottom] = dirichlet(T_BC);
f[bottom] = neumann(0);




event init (t = 0) {
    if (!restore (file = "restart")) {
        int iter = 0;
        do {
            iter++;
            foreach()
            {
                T[] = TMIN;
                f[] = (sq(x+2/8.) + sq(y) - sq(0.25/8.0) > 0 &&
                       sq(x+3.2/8.0) + sq(y-2.0/8.0) - sq(0.25/8.0) > 0 &&
                       sq(x+3/8.0) + sq(y-1.0/8.0) - sq(0.3/8.0) > 0 &&
                       sq(x+2.6/8.0) + sq(y+1.5/8.0) - sq(0.3/8.0) > 0 &&
                       sq(x+2.8/8.0) + sq(y+3./8.0) - sq(0.03*L0) > 0) ? 1 : 0;// && (x < 2)

                fs[] = sq(x-0.25*L0) + sq(y) - sq(0.03*L0) > 0 ? 0 : 1;
                u.x[] = U_BC*(1.0 - fs[]);
            }
            boundary ({f, fs, T, u});
        }while (adapt_wavelet({f, fs, T, u}, (double []){feps, feps, Teps, ueps, ueps},
                              maxlevel = MAXLEVEL, minlevel=MINLEVEL).nf != 0 && iter <= 15);
        fprintf(stderr, "init refinement iter=%d", iter);
        foreach() {
            alpha_doc[] = 0;
        }
        fprintf(stderr, "kappa1=%g kappa2=%g kappa3=%g ", kappa1,kappa2,kappa3);
        foreach_face() kappav.x[] = fm.x[]*kappav(f[], fs[]);
        advection_centered(T, u, u_grad_scalar);
    }
}

event velocity_correction(i++){
    foreach() foreach_dimension() u.x[] *= (1-fs[]);
}
//event adapt_step(i<=5)  DT = 1e-9;//event adapt_step(i<=5)  DT = 1e-9;

//Output
#include "../src_local/output_vtu_foreach.h"
event vtk_file (i += 1)
{
    char subname[80]; sprintf(subname, "hrhs");
    scalar l[]; foreach() l[] = level;
    //be careful with kappa, mu. They can be const unity

//    vector mapped_kappa_lower[], mapped_kappa_upper[];
//    vector mapped_mu_lower[], mapped_mu_upper[];
//    face_vector2vector(kappa, mapped_kappa_lower, mapped_kappa_upper);
//    face_vector2vector(kappa, mapped_mu_lower, mapped_mu_upper);


//    vector mapped_kappa_lower[], mapped_kappa_upper[];
//    vector mapped_mu_lower[], mapped_mu_upper[];
//    foreach()
//    foreach_dimension(){
//        mapped_kappa_lower.x[] = kappa.x[];
//        mapped_kappa_upper.x[] = kappa.x[1];
//        mapped_mu_lower.x[] = mu.x[];
//        mapped_mu_upper.x[] = mu.x[1];
//    }
    advection_centered(T, u, u_grad_scalar);// DELETE THIS LINE!
    output_vtu_MPI( (scalar *) {l, f, fs, T, alpha_doc, rho, u_grad_scalar}, (vector *) {u, kappa, mu}, subname);
}

#if DUMP
event snapshot (i += 1000)
{
  char name[80];
  sprintf (name, "dump-%d", i);
  dump (file = name);
}
#endif

event adapt (i+=100) {
    adapt_wavelet ({f, fs, T, alpha_doc, u}, (double[]){feps, feps, Teps, aeps, ueps, ueps}, MAXLEVEL, MINLEVEL);
}

event stop(t=200);

/**
We check the number of iterations of the Poisson and viscous
problems. */

//event logfile (i++)
//fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);


/**
We produce animations of the vorticity and tracer fields... */

//event movies (t+=0.01; i <= 5000.)
//{
//    scalar omega[];
//    vorticity (u, omega);
//
//    output_ppm (omega, file = "vort.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
//            min = -10, max = 10, linear = true);
//    output_ppm (f, file = "f.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
//            linear = false, min = 0, max = 1);
//}

