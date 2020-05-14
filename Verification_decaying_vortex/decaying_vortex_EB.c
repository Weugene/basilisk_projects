/**
# Decaying vortex problem
To verify the spatial accuracy of the present method, the decaying vortex problem is chosen because it is an unsteady problem with an analytical solution:
 $$u(x, y, t) = − \cos 􏱲x \sin 􏱲y \exp{−2􏱲\pi^2 t/Re}$$
 $$v(x, y, t) =   \sin 􏱲x \cos 􏱲y \exp{−2􏱲\pi^2 t/Re}$$
 $$p(x, y, t) =  −\frac14(\cos 2􏱲x + \cos 2􏱲y)\exp{−4\pi^2 t/Re}$$

The computational domain is −1.5 < x, y 􏰒< 1.5 and the IB is located at x = ± 1 and y = ± 1.
The Reynolds number based on the maximum velocity and vortex size is set to 30, and the initial and boundary conditions are given by the analytical solution above.
Simulations are performed till $t=0.3$.
We use the centered Navier-Stokes solver, with embedded boundaries and
advect the passive tracer *f*. */

#define FILTERED
#define RELATIVE_RESIDUAL
#define EPS_MAXA 2
#include "embed.h"
#include "navier-stokes/centered.h"
// #include "navier-stokes/perfs.h"
#include "../src_local/output_vtu_foreach.h"
#include "../src_local/utils-weugene.h"
face vector muv[];
vector frame_normal[], uexact[];
scalar omega[], utau[];
scalar pexact[];
int maxlevel = 8;
int minlevel = 4;
double Ldomain = 3.0, RE = 30.;

void frame_embed(scalar cs, face vector fs, vector frame_normal, double t){
    vertex scalar phi[];
    foreach_vertex() {
        phi[] = ( fabs(x) <= 1 && fabs(y) <= 1 ) ? 1 : -1;
    }
    boundary ({phi});
    fractions (phi, cs, fs);
    double un;
    foreach() {
        foreach_dimension() frame_normal.x[] = (fs.x[1] - fs.x[])/Delta;
        un = 0;
        foreach_dimension() un += sq(frame_normal.x[]);
        un = sqrt(un + SEPS);
        foreach_dimension() frame_normal.x[] /= un;
    }

}

#define uxe (-cos(pi*x) * sin(pi*y) * exp(-2.0 * sq(pi) * t / RE))
#define uye ( sin(pi*x) * cos(pi*y) * exp(-2.0 * sq(pi) * t / RE))
#define pxe (-0.25 *  (cos(2*pi*x) + cos(2.0*pi*y)) * exp(-4.0 * sq(pi) * t / RE))
void theory(vector u, scalar p, vector n, double RE){
    double un;
    foreach(){
        u.x[] = -cos(pi*x) * sin(pi*y) * exp(-2.0 * sq(pi) * t / RE);
        u.y[] =  sin(pi*x) * cos(pi*y) * exp(-2.0 * sq(pi) * t / RE);
        p[] = -0.25 *  (cos(2*pi*x) + cos(2.0*pi*y)) * exp(-4.0 * sq(pi) * t / RE);
        un = u.x[]*n.x[] + u.y[]*n.y[];
        utau[] = sqrt( sq(u.x[] - un*n.x[]) + sq(u.y[] - un*n.y[]) );
    }
}
//$$u(x, y, t) = − \cos 􏱲x \sin 􏱲y \exp{−2􏱲\pi^2 t/Re}$$
//$$v(x, y, t) =   \sin 􏱲x \cos 􏱲y \exp{−2􏱲\pi^2 t/Re}$$
//$$p(x, y, t) =  −\frac14(\cos 2􏱲x + \cos 2􏱲y)\exp{−4\pi^2 t/Re}$$


#define unor (uxe*frame_normal.x[] + uye*frame_normal.y[])
#define utauval sqrt( sq(uxe - (uxe*frame_normal.x[] + uye*frame_normal.y[])*frame_normal.x[]) + sq(uye - (uxe*frame_normal.x[] + uye*frame_normal.y[])*frame_normal.y[]) )

u.n[left]  = dirichlet(unor);
u.t[left]  = dirichlet(utau[]);
p[left]    = dirichlet(pxe);
pf[left]   = dirichlet(pxe);

u.n[right] = dirichlet(unor);
u.t[right] = dirichlet(utau[]);
p[right]   = dirichlet(pxe);
pf[right]  = dirichlet(pxe);

u.n[top] = dirichlet(unor);
u.t[top] = dirichlet(utau[]);
p[top]   = dirichlet(pxe);
pf[top]  = dirichlet(pxe);

u.n[bottom] = dirichlet(unor);
u.t[bottom] = dirichlet(utau[]);
p[bottom]   = dirichlet(pxe);
pf[bottom]  = dirichlet(pxe);

/**
All walls are no-slip. */

u.n[embed] = dirichlet(unor);
u.t[embed] = dirichlet(utau[]);


int main(int argc, char * argv[]) {
    if (argc > 1) {
        maxlevel = atoi(argv[1]); //convert from string to float
    }
    size (Ldomain);
    origin (-0.5*Ldomain, -0.5*Ldomain);
//    DT = 1e-8;
    DT = 1e-4;
    CFL = 0.4;
    TOLERANCE = 1e-8;
    RELATIVE_RES_TOLERANCE = 0.1;
    NITERMAX = 30;
    N = 512;
    mu = muv;

    run();
}

/**
We set a constant viscosity corresponding to a Reynolds number of 40, 100,
based on the cylinder diameter (1) and the inflow velocity (1). */

event properties (i++)
{
    foreach_face() muv.x[] = fm.x[]/RE;
}


event init (t = 0)
{
    if (!restore (file = "restart")) {
        int it = 0;
        do {
            frame_embed (cs, fs, frame_normal, 0);
            theory(u, p, frame_normal, utau, RE);
            boundary({u.x, u.y});
            foreach() {
                if (norm(frame_normal)>2*sqrt(SEPS)) fprintf(ferr, "ada %g %g \n", frame_normal.x[], frame_normal.y[]);
            }
        } while ( ++it <= 10 && adapt_wavelet({cs, u}, (double []){1e-5, 1e-4, 1e-4}, maxlevel=maxlevel, minlevel=minlevel).nf != 0);
        foreach() {
            double un = 0;
            foreach_dimension() un += sq(frame_normal.x[]);
            un = sqrt(un + SEPS);
            foreach_dimension() frame_normal.x[] /= un;
            if (norm(frame_normal)>2*sqrt(SEPS)) fprintf(ferr, "it=%d w %g %g \n", it, frame_normal.x[], frame_normal.y[]);
        }
        event("vtk_file");
    }

}

event set_dtmax (i++) {
    if (i<=100) {
        NITERMIN=100;
        NITERMAX=150;
    }else{
        NITERMIN=10;
        NITERMAX=30;
    }
    DT *= 1.05;
    DT = min(DT, CFL*Ldomain/pow(2, maxlevel+3));
    fprintf(ferr, "set_dtmax: tnext= %g Dt= %g", tnext, DT);
}

/**
We check the number of iterations of the Poisson and viscous
problems. */

event logfile (i++){
    fprintf (ferr, "i=%d t=%g dt=%g iter_p=%d iter_u=%d \n", i, t, dt, mgp.i, mgu.i);
}
/**
We produce animations of the vorticity and tracer fields... */

//event movies (i += 4; t <= 15.)
//{
//    scalar omega[], m[];
//    vorticity (u, omega);
//    foreach()
//            m[] = cs[] - 0.5;
//    boundary ({m});
//    output_ppm (omega, file = "vort.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
//        min = -10, max = 10, linear = true, mask = m);
//#if TURN_ON_TRACER == 1
//    output_ppm (f, file = "f.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
//        linear = false, min = 0, max = 1, mask = m);
//#endif
//}

//Output
event vtk_file (i++){
//event vtk_file (t += 1){
    char subname[80]; sprintf(subname, "vortex_EB");
    scalar l[], omega[];
    vorticity (u, omega);
    foreach() {l[] = level;}
    theory(uexact, pexact, frame_normal, utau, RE);

    output_vtu_MPI( (scalar *) {cs, omega, p, pexact, l}, (vector *) {u, uexact, frame_normal}, subname, 0 );
}
/**
We adapt according to the error on the embedded geometry, velocity*/
#define ADAPT_SCALARS {cs, u}
#define ADAPT_EPS_SCALARS {1e-5, 1e-3, 1e-3}
event adapt (i++){
    double eps_arr[] = ADAPT_EPS_SCALARS;
    MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
}

event stop (t = 0.3);
/**
## See also

* [Same example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/cylinder.html)
*/
