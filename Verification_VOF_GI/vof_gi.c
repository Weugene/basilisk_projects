
//#define FILTERED
//#define JACOBI 1

scalar omega[];
#include "navier-stokes/centered.h"
#include "view.h"
#include "../src_local/output_vtu_foreach.h"
#include "two-phase.h"
#include "tension.h"

int maxlevel = 10;
int minlevel = 4;
double U0=1, rhol=1, sig=0.0005, Lchar=1, mul=1, Lb=0.3, xs0 = -0.4, Rb=0.0625, rad=0.0625;
double Rrho=1, Rmu=1;
double RE, CA, Rrad;


/**
The domain is the periodic unit square centered on the origin. */

/**
The fluid is injected on the left boundary with a unit velocity. The
tracer is injected in the lower-half of the left boundary. An outflow
condition is used on the right boundary. */

//u.n[left]  = dirichlet(0);
//u.t[left]  = dirichlet(0);
//p[left]    = neumann(0);
//pf[left]   = neumann(0);
//
//u.n[right] = neumann(0.);
//u.t[right] = neumann(0.);
//p[right]   = neumann(0.);
//pf[right]  = neumann(0.);

int main(int argc, char * argv[])
{
    RE = U0*Lchar*rhol/mul;
    CA = U0*mul/sig;
    Rrad = Rb/rad;
    if (argc > 1) {
        RE = atof(argv[1]);
    }
    if (argc > 2) {
        CA = atof(argv[2]);
    }
    if (argc > 3) {
        Rrad = atof(argv[3]);
    }
    if (argc > 4) {
        Rrho = atof(argv[4]);
    }
    if (argc > 5) {
        Rmu = atof(argv[5]);
    }
    if (argc > 6) {
        maxlevel = atoi(argv[6]);
    }
    size (1.0);
    origin (-0.5*L0, -0.5*L0);

    TOLERANCE = 1e-4;
    NITERMAX = 30;
    N = 1 << minlevel;
    periodic(top);
    periodic(right);

    Rb = Rrad*0.0625;
    rho1 = 1.; rho2 = rho1/Rrho;
    mu1 = 1./RE; mu2 = mu1/Rmu;
    f.sigma = 1./RE/CA;
    fprintf(ferr, "maxlevel=%d tol=%g NITERMAX=%d\n"
                  "RE=%g CA=%g Rb/rad=%g rho1/rho2=%g mu1/mu2=%g\n"
                  "mu1=%g mu2=%g rho1=%g rho2=%g sigma=%g\n",
            maxlevel, TOLERANCE, NITERMAX,
            RE, CA, Rrad, Rrho, Rmu,
            mu1, mu2, rho1, rho2, f.sigma);


    run();
}

scalar divu[];
void bubbles (scalar f)
{
    const int ns=1;
    double xc[ns], yc[ns], R[ns];
    xc[0] = Lb + xs0; yc[0] = 0; R[0] = Rb;

    vertex scalar phi[];
    face vector ff[];
    foreach_vertex() {
        phi[] = HUGE;
        for (double xp = -L0; xp <= L0; xp += L0)
            for (double yp = -L0; yp <= L0; yp += L0)
                for (int i = 0; i < ns; i++)
                    phi[] = intersection (phi[], (sq(x + xp - xc[i]) + sq(y + yp - yc[i]) - sq(R[i])));
        //phi[] = -phi[];
    }
    
    boundary ({phi});
    fractions (phi, f, ff);
}
event init (t = 0) {
    if (!restore (file = "restart")) {
        int it = 0;
        do {
            bubbles(f);
            foreach() {
                u.x[] = 1;
                u.y[] = 0;
            }
            boundary({f,u});
        }while (adapt_wavelet({f, u}, (double []){1e-5, 1e-3, 1e-3}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
        DT=1e-10;
    }
}

event set_dtmax (i++) if (i<500) DT *= 1.05;

event vof(i++){
    foreach() { u.x[]=1; u.y[]=0; }
    foreach_face(x) uf.x[] = 1;
    foreach_face(y) uf.y[] = 1e-9;
}

void correct_press(scalar p, int i){
    double press = 0;
    int ip = 0;
#if 1 // Left bottom Corner
    foreach(){
        if (ip == 0){
            press = p[];
            ip++;
            @if _MPI
                MPI_Bcast(&press, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            @endif
        }
    }
#else //average value
    press = normf(p).avg;
#endif
    foreach(){
        p[] -= press;
    }
}

event end_timestep(i++){
    correct_press(p, i);
    vorticity (u, omega);
}

event logfile (i++) {
    double avggas = sq(L0) - normf(f).avg;
    foreach() {
        divu[] = 0;
        foreach_dimension() divu[] += (uf.x[1]-uf.x[])/Delta;
    }
    double Linf_u = -10, Linf_omega_min = 1e10, Linf_omega_max = -10;
    foreach( reduction(max:Linf_u) reduction(min:Linf_omega_min) reduction(max:Linf_omega_max)){
        if (fabs(divu[]) > Linf_u) Linf_u = fabs(divu[]);
        if (fabs(omega[]) < Linf_omega_min) Linf_omega_min = fabs(omega[]);
        if (fabs(omega[]) > Linf_omega_max) Linf_omega_max = fabs(omega[]);
    }
    fprintf (ferr, "i=%d t=%g dt=%g iter_p=%d iter_u=%d AvgGas=%g divu=%g Omega_min=%g Omega_max=%g\n", i, t, dt, mgp.i, mgu.i, avggas, Linf_u, Linf_omega_min, Linf_omega_max);
}

/**
We produce animations of the vorticity and tracer fields... */

event images (t += 0.1) {
    static FILE * fp = popen ("ppm2gif > vort.gif", "w");
    output_ppm (omega, fp, min=-10, max=10, linear=true);
}

//Output
//event vtk_file (i>1000000){
//event vtk_file (i += 1){
event vtk_file (t += 0.01){
    char subname[80]; sprintf(subname, "myvof");
    scalar l[];
    foreach() {l[] = level;}

    vector mapped_data_lower[], mapped_data_upper[], mu_low[];
    foreach() {
        foreach_dimension()
        {
            mapped_data_lower.x[] = uf.x[];
            mapped_data_upper.x[] = uf.x[1];
            mu_low.x[]=mu.x[];
        }
    }
    output_vtu_MPI( (scalar *) {f, omega, p, l, divu, rho}, (vector *) {u, g, a, mu_low, mapped_data_lower, mapped_data_upper}, subname, 1 );
}

#define ADAPT_SCALARS {f}
#define ADAPT_EPS_SCALARS {1e-5}
event adapt (i++){
    double eps_arr[] = ADAPT_EPS_SCALARS;
//    MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);

}

event snapshot (t += 0.1) {
      char name[80];
      sprintf (name, "snapshot-%g", t);
      dump (name);
}

event stop(t = 10);