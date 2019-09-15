/**
# Bénard–von Kármán Vortex Street for flow around a cylinder at Re=160

An example of 2D viscous flow around a simple solid boundary. Fluid is
injected to the left of a channel bounded by solid walls with a slip
boundary condition. A passive tracer is injected in the bottom half of
the inlet.

![Animation of the vorticity field.](karman/vort.mp4)(loop)

![Animation of the tracer field.](karman/f.mp4)(loop)

We use the centered Navier-Stokes solver and
advect the passive tracer *f*. */


#include "navier-stokes/centered.h"
//#include "../src_local/three-phase-rheology.h"
#include "two-phase.h"
double rho3 = 1., mu3 = 0.;
scalar alpha_doc[];
scalar T[];
scalar fs[];
// #include "navier-stokes/perfs.h"
#include "tracer.h"
//#include "diffusion.h"


scalar * tracers = {f};

/**
The domain is eight units long, centered vertically. */

int main() {
    L0 = 8.;
    origin (-0.5, -L0/2.);
    N = 1024;
//    mu1=0.1; mu2=5.419e-5; mu3=0.1; // kg/(m*s)
//    rho1=1140; rho2=1.14; rho3=1900; // kg/m3
    mu1=1; mu2=1; mu3=1; // kg/(m*s)
    rho1=1; rho2=1; rho3=1; // kg/m3
    run();
}

/**
The fluid is injected on the left boundary with a unit velocity. The
tracer is injected in the lower-half of the left boundary. An outflow
condition is used on the right boundary. */

u.n[left]  = dirichlet(1);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(y < 0);//a lower half
//muv[left]  = dirichlet((y < 0)?mu1:mu2);//a lower half

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
The top and bottom walls are free-slip and the obstacles are no-slip. */

event init (t = 0) {
    foreach()
    {
        alpha_doc[] = 0;
        T[] = 1;
        u.x[] = 1;
        f[] = 0;
        fs[] = 0;
    }
}

double Htr = 0;
double Arrhenius_const = 4.53e+7;//1/s
double Ea_by_R = 72900/8.314;// Kelvin
double n_degree = 1.667;
double m_degree = 0.333;


//double Cp1 = 1100, Cp2 = 1006, Cp3 = 840;//J/(kg*K)
//double kappa1 = 1100, kappa2 = 0.02535, kappa3 = 1.11;//W/(m*K)
double Cp1 = 1, Cp2 = 1, Cp3 = 1;//J/(kg*K)
double kappa1 = 1, kappa2 = 0.1, kappa3 = 1;//W/(m*K)
void average(scalar A, scalar f, scalar fs, double A1, double A2, double A3){
    foreach(){
        A[] = f[]*(A1 - A2) + A2 + fs[]*(A3 - A2);
    }
}
void advection_centered (scalar f, vector u, scalar df)
{
  foreach()
    df[] = ((f[] + f[-1,0])*u.x[] -
	    (f[] + f[1,0])*u.x[1,0] +
	    (f[] + f[0,-1])*u.y[] -
	    (f[] + f[0,1])*u.y[0,1])/(2.*Delta);
  boundary ((scalar*) {df});
}

void advection_upwind (scalar f, vector u, scalar df)
{
  foreach()
    df[] = ((u.x[] < 0. ? f[] : f[-1,0])*u.x[] -
	    (u.x[1,0] > 0. ? f[] : f[1,0])*u.x[1,0] +
	    (u.y[] < 0. ? f[] : f[0,-1])*u.y[] -
	    (u.y[0,1] > 0. ? f[] : f[0,1])*u.y[0,1])/Delta;
    boundary ((scalar*) {df});
}


mgstats mgT;
//event advance_alpha_T (i++){
//    //dt = dtnext (0.01);
//    scalar r[], theta[];
//    scalar Cp[];
//    face vector kappa[];
//    average(theta, f, fs, rho1*Cp1, rho2*Cp2, rho3*Cp3);
//    foreach(){
//        kappa.x[] = f[]*(kappa1 - kappa2) + kappa2 + fs[]*(kappa3 - kappa2);
//    }
//
//    scalar u_grad_scalar[], tmp[];
//    advection_centered(T, u, u_grad_scalar);
//    foreach() {
//        tmp[] = Arrhenius_const*pow(1-alpha_doc[], n_degree)*exp(-Ea_by_R/T[]);
//        r[] = Htr*rho1*f[]*tmp[] - theta[]*u_grad_scalar[];
//    }
//    mgT = diffusion (T, dt, D = kappa, r = r, theta = theta);
//
//    fprintf (stderr, "mg: i=%d t=%g p=%d u=%d T=%d\n", i, t, mgp.i, mgu.i, mgT.i);
//    advection_centered(alpha_doc, u, u_grad_scalar);
//
//    foreach() {
//        alpha_doc[] += dt*(tmp[] - u_grad_scalar[]);
//    }
//}

//Output
static int iteration=0;
#include "output_fields/output_vtu_foreach.h"
event vtk_file (i++)
{
    int nf = iteration;
    scalar l[];
    foreach()
        l[] = level;

    char name[80], subname[80];
    FILE *fp;
    sprintf(name, "hrhs_%4.4d_n%3.3d.vtu", nf, pid());
    fp = fopen(name, "w");

    output_vtu_bin_foreach((scalar *) {l, f, fs, rho, T, alpha_doc}, (vector *) {u, mu}, 64, fp, false);
    fclose(fp);
    @if _MPI
        if (pid() == 0) {
            sprintf(name, "hrhs_%4.4d.pvtu", nf);
            sprintf(subname, "hrhs_%4.4d", nf);
            fp = fopen(name, "w");
            output_pvtu_bin((scalar *) {l, f, fs, rho, T, alpha_doc}, (vector *) {u, mu}, 64, fp, subname);
            fclose(fp);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    @endif
    fprintf (ferr, "iteration: %d\n", iteration); fflush (ferr);
    iteration++;
}

#if DUMP
event snapshot (i += 1000)
{
  char name[80];
  sprintf (name, "dump-%d", i);
  dump (file = name);
}
#endif

/**
We adapt according to the error on the embedded geometry, velocity and
tracer fields. */

event adapt (i++) {
    adapt_wavelet ({f}, (double[]){1e-3}, 11, 4);
}

event stop(i = 50);
/**
## See also

* [Same example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/cylinder.html)
*/





/**
We check the number of iterations of the Poisson and viscous
problems. */

//event logfile (i++)
//fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

//face vector av[];
//event acceleration (i++) {
//    foreach_face()
//        av.x[] = (cs[]/(0.01))*((u.x[] + u.x[-1])/2. );
//    foreach_face()
//        av.x[] -= cs[]*(u.x[]);
//    boundary ((scalar *){av});
//}

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
































































/**
# Bénard–von Kármán Vortex Street for flow around a cylinder at Re=160

An example of 2D viscous flow around a simple solid boundary. Fluid is
injected to the left of a channel bounded by solid walls with a slip
boundary condition. A passive tracer is injected in the bottom half of
the inlet.

![Animation of the vorticity field.](karman/vort.mp4)(loop)

![Animation of the tracer field.](karman/f.mp4)(loop)

We use the centered Navier-Stokes solver, with embedded boundaries and
advect the passive tracer *f*. */


#include "navier-stokes/centered.h"
// #include "navier-stokes/perfs.h"
#include "two-phase.h"
//#include "tracer.h" //don't use it if you added two-phase
#include "diffusion.h"


#define MAXLEVEL 9
#define MINLEVEL 4
#define feps 1e-2
#define Teps 1e-2
#define aeps 1e-2
#define ueps 1e-2
//scalar f[];
scalar fs[]; //comment for 3-phase.h
scalar T[]; //comment for 3-phase.h
scalar alpha_doc[]; //comment for 3-phase.h
face vector kappa[];
double rho3; //comment for 3-phase.h
//scalar * tracers = {f}; //?
//face vector muv[];

int main() {
    L0 = 8.;
    origin (-L0/2, -L0/2.);
    N = 64;
//    mu = muv;
    run();
    rho1=1; rho2=0.1; rho3=1;
    mu1=1; mu2=0.1;

}
u.n[left]  = dirichlet(1);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(1);
T[left]    = dirichlet(1);
fs[left]   = dirichlet(0);
alpha_doc[left] = dirichlet(0);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);
T[right]    = neumann(0);
fs[right]   = neumann(0);
alpha_doc[right] = neumann(0);

#define UT_BC exp(-pow(x+L0/2.,2)/(2*pow(L0/10.,2)))
#define T_BC 1.5+0.5*tanh(x/(L0/10))
u.n[top] = dirichlet(0);
u.t[top] = neumann(0);
//u.t[top] = dirichlet(UT_BC);
alpha_doc[top] = neumann(0);
T[top] = dirichlet(T_BC);
fs[top] = neumann(0);
f[top] = neumann(0);

u.n[bottom] = dirichlet(0.);
u.t[bottom] = neumann(0);
//u.t[bottom] = dirichlet(UT_BC);
alpha_doc[bottom] = neumann(0);
T[bottom] = dirichlet(T_BC);
f[bottom] = neumann(0);

/**
The top and bottom walls are free-slip and the cylinder is no-slip. */
double Htr = 0;
//double Arrhenius_const = 4.53e+7;//1/s
//double Ea_by_R = 72900/8.314;// Kelvin
double Arrhenius_const = 10;//1/s
double Ea_by_R = 10;// Kelvin
double n_degree = 1.667;
double m_degree = 0.333;


//double Cp1 = 1100, Cp2 = 1006, Cp3 = 840;//J/(kg*K)
//double kappa1 = 1100, kappa2 = 0.02535, kappa3 = 1.11;//W/(m*K)
double Cp1 = 1, Cp2 = 1, Cp3 = 1;//J/(kg*K)
double kappa1 = 1, kappa2 = 0.1, kappa3 = 1;//W/(m*K)

event init (t = 0) {
    if (!restore (file = "restart")) {
        int iter = 0;
        do {
            iter++;
            foreach()
            {
                T[] = 1;
                f[] = (sq(x) + sq(y) - sq(1) > 0) ? 1 : 0;
                u.x[] = 1;
            }
            boundary ({f, T, u});
        }while (adapt_wavelet({f, T, u}, (double []){feps, Teps, ueps, ueps},
                              maxlevel = MAXLEVEL, minlevel=MINLEVEL).nf != 0 && iter <= 10);
        foreach() {
            alpha_doc[] = 0;
            fs[] = 0;
        }
        foreach_face(){
            kappa.x[] = f[]*(kappa1 - kappa2) + kappa2 + fs[]*(kappa3 - kappa2);
        }
    }
}

event adapt_step(i<=5)  DT = 1e-9;//event adapt_step(i<=5)  DT = 1e-9;


void average(scalar A, const scalar f, const scalar fs, const double A1, const double A2, const double A3){
    foreach(){
        A[] = f[]*(A1 - A2) + A2 + fs[]*(A3 - A2);
    }
}
void advection_centered (scalar f, vector u, scalar df)
{
    foreach()
    df[] = ((f[] + f[-1,0])*u.x[] -
            (f[] + f[1,0])*u.x[1,0] +
            (f[] + f[0,-1])*u.y[] -
            (f[] + f[0,1])*u.y[0,1])/(2.*Delta);
    boundary ((scalar*) {df});
}

void advection_upwind (scalar f, vector u, scalar df)
{
    foreach()
    df[] = ((u.x[] < 0. ? f[] : f[-1,0])*u.x[] -
            (u.x[1,0] > 0. ? f[] : f[1,0])*u.x[1,0] +
            (u.y[] < 0. ? f[] : f[0,-1])*u.y[] -
            (u.y[0,1] > 0. ? f[] : f[0,1])*u.y[0,1])/Delta;
//    boundary ((scalar*) {df});
}


mgstats mgT;
event advance_alpha_T (i++,last){
dt = dtnext (0.01);//??
scalar r[], theta[];

average(theta, f, fs, rho1*Cp1, rho2*Cp2, rho3*Cp3);
foreach_face(){
    kappa.x[] = f[]*(kappa1 - kappa2) + kappa2 + fs[]*(kappa3 - kappa2);
}

scalar u_grad_scalar[], tmp[];
advection_centered(T, u, u_grad_scalar);
//    advection_upwind (T, u, u_grad_scalar);
foreach() {
    tmp[] = Arrhenius_const*pow(1-alpha_doc[], n_degree)*exp(-Ea_by_R/T[]);
    r[] = Htr*rho1*f[]*tmp[] - theta[]*u_grad_scalar[];
}
mgT = diffusion (T, dt, D = kappa, r = r, theta = theta);

fprintf (stderr, "mg: i=%d t=%g p=%d u=%d T=%d\n", i, t, mgp.i, mgu.i, mgT.i); //number of iterations
advection_centered(alpha_doc, u, u_grad_scalar);
//    advection_upwind (alpha_doc, u, u_grad_scalar);
foreach() {
    alpha_doc[] += dt*(tmp[] - u_grad_scalar[]);
}
boundary ((scalar*) {alpha_doc});

}
/**
We check the number of iterations of the Poisson and viscous
problems. */

//event logfile (i++)
//fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

//event acceleration (i++) {
//    foreach_face()
//        av.x[] = (cs[]/(0.01))*((u.x[] + u.x[-1])/2. );
//    foreach_face()
//        av.x[] -= cs[]*(u.x[]);
//    boundary ((scalar *){av});
//}

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

//Output
static int iteration=0;
#define OUTPUT_VARS (scalar *) {l, f, T, alpha_doc}, (vector *) {u, mu}
#include "output_fields/output_vtu_foreach.h"
event vtk_file (i++)
{
int nf = iteration;
scalar l[];
foreach()
        l[] = level;
        char name[80], subname[80];
        FILE *fp;
                sprintf(name, "hrhs_%4.4d_n%3.3d.vtu", nf, pid());
fp = fopen(name, "w");
output_vtu_bin_foreach(OUTPUT_VARS, 64, fp, false);
fclose(fp);

if (pid() == 0) {
sprintf(name, "hrhs_%4.4d.pvtu", nf);
sprintf(subname, "hrhs_%4.4d", nf);
fp = fopen(name, "w");
output_pvtu_bin(OUTPUT_VARS, 64, fp, subname);
fclose(fp);
}
@if _MPI
MPI_Barrier(MPI_COMM_WORLD);
@endif
fprintf (ferr, "iteration: %d\n", iteration); fflush (ferr);
iteration++;
}




//#if DUMP
//event snapshot (i += 1000)
//{
//  char name[80];
//  sprintf (name, "dump-%d", i);
//  dump (file = name);
//}
//#endif

/**
We adapt according to the error on the embedded geometry, velocity and
tracer fields. */

event adapt (i++) {
adapt_wavelet ({f, T, alpha_doc, u}, (double[]){feps, Teps, aeps, ueps, ueps}, MAXLEVEL, MINLEVEL);
}

event stop(i = 10000);
/**
## See also

* [Same example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/cylinder.html)
*/