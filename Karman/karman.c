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


#include "../src_local/centered-weugene.h"
// #include "navier-stokes/perfs.h"
//#include "two-phase.h"
#include "../src_local/three-phase-rheology.h"
#include "tension.h"
//#include "tracer.h" //don't use it if you added two-phase
#include "diffusion.h"



//#define TMAX 400
//#define TMIN 300

#define MAXLEVEL 10
#define MINLEVEL 4
#define feps 1e-3
#define Teps 1e-3
#define aeps 1e-3
#define ueps 1e-3
#undef SEPS
#define SEPS 1e-10

#define Re 1e+0//,lmlm/
#define Pe 1e+0
#define Ec 1e+0
#define Ex 1e+0
#define Ar 5e+0
#define Rrho 1e+2
#define Rmu 1e+2
#define Rkappa 1e+2
#define RT 2e+0

#define SIGMA 1e+0
#define MAXVEL 0.2
#define TMAX 2
#define TMIN 1

face vector kappa[];

int main() {
    L0 = 8.;
    origin (-L0/2, -L0/2.);
    N = 512;
    CFL = 0.1;
    DT = 1e-4;
//    stokes = true;
//    rho1 = 1140; rho2 = 1; rho3 = 2000;
//    mu1 = 0.155; mu2 = 1.81e-5; mu3 = 1;
    rho1 = 1.0/Re; rho2 = 1.0/(Re*Rrho); rho3 = 1;
    mu1 = 1;  mu2 = 0.1;  mu3 = 1;
    f.sigma = SIGMA; chi = 3;
    run();
}
#define U_BC MAXVEL*(sqrt(1 - pow(y/(L0/2.),16)))
//#define UT_BC exp(-pow(x + L0/2.,2)/(2*pow(L0/10.,2)))
#define T_BC (0.5*(TMAX + TMIN) + 0.5*(TMAX - TMIN)*tanh((x+0.45*L0)/(L0/10)))
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
alpha_doc[right] = neumann_homogeneous();


u.n[top] = dirichlet(0);
u.t[top] = neumann(0);
//u.t[top] = dirichlet(0);
alpha_doc[top] = neumann(0);
T[top] = dirichlet(T_BC);
fs[top] = neumann(0);
f[top] = neumann(0);

u.n[bottom] = dirichlet(0.);
u.t[bottom] = neumann(0);
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
                T[] = 1;
                f[] = (sq(x+3) + sq(y) - sq(0.25) > 0 &&
                       sq(x+3.2) + sq(y-2) - sq(0.25) > 0 &&
                       sq(x+3) + sq(y-1) - sq(0.3) > 0 &&
                       sq(x+2.6) + sq(y+1.5) - sq(0.3) > 0 &&
                       sq(x+2.8) + sq(y+3) - sq(0.4) > 0) && (x < -2) ? 1 : 0;
                u.x[] = U_BC;
            }
            boundary ({f, T, u});
        }while (iter <=5 || adapt_wavelet({f, T, u}, (double []){feps, Teps, ueps, ueps},
                              maxlevel = MAXLEVEL, minlevel=MINLEVEL).nf != 0 && iter <= 15);
        fprintf(stderr, "init refinement iter=%d", iter);
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

//Output
event vtk_file (t += 0.01)
{
    char subname[80]; sprintf(subname, "hrhs");
    //be careful with kappa, mu. They can be const unity
    output_vtu_MPI( (scalar *) {l, f, T, alpha_doc, thetav, r, rho, u_grad_scalar, tmp}, (vector *) {u, kappa, mu}, subname)
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
    adapt_wavelet ({f, T, alpha_doc, u}, (double[]){feps, Teps, aeps, ueps, ueps}, MAXLEVEL, MINLEVEL);
//    unrefine (x > 0.45*L0 && x<-0.45*L0);

}

event stop(t=100);
/**
## See also

* [Same example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/cylinder.html)
*/