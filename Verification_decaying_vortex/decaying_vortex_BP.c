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

#define BRINKMAN_PENALIZATION 1
#define DEBUG_BRINKMAN_PENALIZATION 1
//#define DEBUG_OUTPUT_VTU_MPI
#define DEBUG_MINMAXVALUES
#define DEBUG_MODE
#define FILTERED
#define JACOBI 1
#define EPS_MAXA 2
//#define RELATIVE_RESIDUAL
#define PERIODIC_BC
scalar omega[], fs[];
vector target_Uv[];
vector deltag[];
scalar deltap[];
scalar phi[];
int maxlevel;
int minlevel;
//int cells_per_zone;
double Ldomain = 4.0, RE = 30., rho1 = 1, mu1;
double eps = 1e-2;
//double dx_min;

#ifdef DEBUG_MODE
    scalar divutmp[], divutmpAfter[];
    scalar mod_du_dx[];
    vector conv_term[], tmp_err_u[];
#endif

#ifdef PERIODIC_BC
    int periodic_BC = 1;
#else
    int periodic_BC = 0;
#endif
#include "../src_local/centered-weugene.h"
//#include "view.h"
#include "../src_local/output_vtu_foreach.h"
#include "fractions.h"
face vector muv[];

void frame_BP(scalar fs, double t){
    vertex scalar phi[];
    face vector fs_face[];
    foreach_vertex() {
        phi[] = ( fabs(x) <= 1 && fabs(y) <= 1 ) ? 0 : 1;
    }
    boundary ({phi});
    fractions (phi, fs, fs_face);
}

#define uxe(x,y,t) (-cos(pi*(x)) * sin(pi*(y)) * exp(-2.0 * sq(pi) * (t) / RE))
#define uye(x,y,t) ( sin(pi*(x)) * cos(pi*(y)) * exp(-2.0 * sq(pi) * (t) / RE))
#define pe(x,y,t)  (-0.25 * (cos(2*pi*(x)) + cos(2.0*pi*(y))) * exp(-4.0 * sq(pi) * (t) / RE))
#define u0 (exp(-2.0 * sq(pi) * (t) / RE))
#define p0 (-0.5 * exp(-4.0 * sq(pi) * (t) / RE))

//#define dpdx(x_,y_,t_,shift_) 0
//#define dpdy(x_,y_,t_,shift_) 0
//#define dpdx(x,y,t) (0.5 * pi * sin(2*pi*x) * exp(-4.0 * sq(pi) * t / RE))
//#define dpdy(x,y,t) (0.5 * pi * sin(2*pi*y) * exp(-4.0 * sq(pi) * t / RE))
//#define dpdx(x,y,t) rho1*( -(uxe(x,y,t) - uf.x[ghost])/dt + (uf.x[ghost] - uxe(x,y,t))/eta_s)
//#define dpdy(x,y,t) rho1*( -(uye(x,y,t) - uf.y[ghost])/dt + (uf.y[ghost] - uye(x,y,t))/eta_s)
//#define dpdx(x_,y_,dt_,shift_) rho1*( -(uxe((x_),(y_),t+dt_) - uf.x[shift_,0])/dt_ )
//#define dpdy(x_,y_,dt_,shift_) rho1*( -(uye((x_),(y_),t+dt_) - uf.y[0,shift_])/dt_ )
#define dpdx(x_,y_,dt_,shift_) rho1*( -(uxe((x_),(y_),t+dt_) - uxe((x_),(y_),t))/dt_ )
#define dpdy(x_,y_,dt_,shift_) rho1*( -(uye((x_),(y_),t+dt_) - uye((x_),(y_),t))/dt_ )

#define dudxe (pi * sin(pi*(x)) * sin(pi*(y)) * exp(-2.0 * sq(pi) * (t) / RE))
#define dudye (pi * sin(pi*(x)) * sin(pi*(y)) * exp(-2.0 * sq(pi) * (t) / RE))

void theory(vector u, scalar p, double t, double RE){
    foreach(){
        u.x[] = uxe(x,y,t);
        u.y[] = uye(x,y,t);
        p[]   = pe(x,y,t);
    }
    boundary({u, p});
    fprintf(ferr, "theory  t= %g\n", t);
}
//$$u(x, y, t) = − \cos 􏱲x \sin 􏱲y \exp{−2􏱲\pi^2 t/Re}$$
//$$v(x, y, t) =   \sin 􏱲x \cos 􏱲y \exp{−2􏱲\pi^2 t/Re}$$
//$$p(x, y, t) =  −\frac14(\cos 2􏱲x + \cos 2􏱲y)\exp{−4\pi^2 t/Re}$$
//in 2D each cells must have 2 BC, in 3D - 3BC
u.n[left]   = dirichlet(uxe(x,y,t));
u.n[right]  = dirichlet(uxe(x,y,t));
u.n[bottom] = dirichlet(uye(x,y,t));
u.n[top]    = dirichlet(uye(x,y,t));

u.t[left]   = dirichlet(uye(x,y,t)); // available x,y,z at the computational domain
u.t[right]  = dirichlet(uye(x,y,t));
u.t[bottom] = dirichlet(uxe(x,y,t));
u.t[top]    = dirichlet(uxe(x,y,t));

//p[left]   = neumann(0);
//p[right]  = neumann(0);
//p[bottom] = neumann(0);
//p[top]    = neumann(0);

p[left]   = dirichlet(pe(x,y,t));
p[right]  = dirichlet(pe(x,y,t));
p[bottom] = dirichlet(pe(x,y,t));
p[top]    = dirichlet(pe(x,y,t));

//pf[left]   = neumann(dpdx(x,y,0.5*dt,0));
//pf[right]  = neumann(dpdx(x,y,0.5*dt,1));
//pf[bottom] = neumann(dpdy(x,y,0.5*dt,0));
//pf[top]    = neumann(dpdy(x,y,0.5*dt,1));

//uf.t[left]   = uye(x,y,t);
//uf.t[right]  = uye(x,y,t);
//uf.t[bottom] = uxe(x,y,t);
//uf.t[top]    = uxe(x,y,t);
//
//uf.n[left]   = uxe(x,y,t);
//uf.n[right]  = uxe(x,y,t);
//uf.n[bottom] = uye(x,y,t);
//uf.n[top]    = uye(x,y,t);

deltap[left]   = neumann(0);
deltap[right]  = neumann(0);
deltap[bottom] = neumann(0);
deltap[top]    = neumann(0);
//deltap[left]   = neumann((uf.x[]      - uxe(x,y,t+dt))/dt);
//deltap[right]  = neumann((uf.x[ghost] - uxe(x,y,t+dt))/dt);
//deltap[bottom] = neumann((uf.y[]      - uye(x,y,t+dt))/dt);
//deltap[top]    = neumann((uf.y[ghost] - uye(x,y,t+dt))/dt);

int main(int argc, char * argv[]) {
    minlevel = 4;
    maxlevel = 13;
    eta_s=1e-8;
//    cells_per_zone = 2;
    eps = 1e-2;
    if (argc > 1) {
        maxlevel = atoi(argv[1]);
    }
    if (argc > 2) {
        eta_s = atof(argv[2]);
    }
    if (argc > 3) {
        eps = atof(argv[3]);
    }
//    dx_min = Ldomain/pow(2,maxlevel);
#ifdef PERIODIC_BC
    periodic(right);
    periodic(top);
#endif
//    eta_s = sq(cells_per_zone*dx_min);
    size (Ldomain);
    origin (-0.5*Ldomain, -0.5*Ldomain);
    DT = 1e-5;
    CFL = 0.3;
    TOLERANCE = 1e-8;
    RELATIVE_RES_TOLERANCE = 0.1;
    NITERMAX = 50;
    mu = muv;
    mu1 = 1.0/RE;
//    reference_pressure = true;
    target_U = target_Uv;
    N = 1<<minlevel;
    fprintf(ferr, "maxlevel= %d eta= %g eps= %g\n", maxlevel, eta_s, eps);
    fprintf(ferr, "TOL=%g NITERMAX=%d Re=%g Rel_res_tol=%g\n", TOLERANCE, NITERMAX, RE, RELATIVE_RES_TOLERANCE);
    fprintf(ferr, "CFL=%g eps=%g\n", CFL, eps);
//    for (int kk = 1; kk <= 4 ; kk++){
//        eta_s *= 10;
        run();
//    }
}

void update_targetU(vector target_Uv, double t){
    foreach() {
        target_Uv.x[] = uxe(x,y,t);
        target_Uv.y[] = uye(x,y,t);
    }
    boundary((scalar *){target_Uv});
}

event init (t = 0)
{
    vector uexact[];
    scalar pexact[];
    if (!restore (file = "restart")) {
        int it = 0;
        do {
            frame_BP(fs, 0);
            theory(uexact, pexact, 0, RE);
        } while (++it <= 10 && adapt_wavelet((scalar *){fs, uexact}, (double []){eps, eps, eps}, maxlevel=maxlevel, minlevel=minlevel).nf != 0);

        foreach() {
            p[] = pexact[];// my initial guess
            pf[] = pexact[];// my initial guess
            foreach_dimension() {
                u.x[] = uexact.x[];
            }
        }
        boundary({p, pf, u});
        update_targetU(target_Uv, 0);
        foreach_face() muv.x[] = mu1;
        boundary((scalar *){muv});
    }

    u.n[left]   = dirichlet(uxe(x,y,t+dt));
    u.n[right]  = dirichlet(uxe(x,y,t+dt));
    u.n[bottom] = dirichlet(uye(x,y,t+dt));
    u.n[top]    = dirichlet(uye(x,y,t+dt));

    u.t[left]   = dirichlet(uye(x,y,t+dt));
    u.t[right]  = dirichlet(uye(x,y,t+dt));
    u.t[bottom] = dirichlet(uxe(x,y,t+dt));
    u.t[top]    = dirichlet(uxe(x,y,t+dt));
}

/**
We set a constant viscosity corresponding to a Reynolds number of 40, 100,
based on the cylinder diameter (1) and the inflow velocity (1). */

event properties (i++) {
    frame_BP(fs, t);
    update_targetU(target_Uv, t+dt);
}

event set_dtmax (i++) {
    if (i==0) event("vtk_file");
    NITERMIN=1;
    NITERMAX=100;

    DT *= 1.05;
//    DT = min(DT, 0.001);
//    DT = min(DT, CFL*Ldomain/pow(2, maxlevel+3));
    fprintf(ferr, "set_dtmax: tnext= %g Dt= %g", tnext, DT);
}

//event advection_term(i++){
//    uf.n[left]   = uxe(x,y,t+0.5*dt);
//    uf.n[right]  = uxe(x,y,t+0.5*dt);
//    uf.n[bottom] = uye(x,y,t+0.5*dt);
//    uf.n[top]    = uye(x,y,t+0.5*dt);
//
//    uf.t[left]   = uye(x,y,t+0.5*dt); // face location at x
//    uf.t[right]  = uye(x,y,t+0.5*dt);
//    uf.t[bottom] = uxe(x,y,t+0.5*dt);
//    uf.t[top]    = uxe(x,y,t+0.5*dt);
//
//    mydt = 0.5*dt;
//    fprintf(ferr, "dt=%g", dt);
//};
//
//event acceleration(i++){
//
//    uf.n[left]   = uxe(x,y,t+dt);
//    uf.n[right]  = uxe(x,y,t+dt);
//    uf.n[bottom] = uye(x,y,t+dt);
//    uf.n[top]    = uye(x,y,t+dt);
//
//    uf.t[left]   = uye(x,y,t+dt); // x, y, is face location
//    uf.t[right]  = uye(x,y,t+dt);
//    uf.t[bottom] = uxe(x,y,t+dt);
//    uf.t[top]    = uxe(x,y,t+dt);
//
//    mydt = dt;
//}

event end_timestep (i += 100){
    vector uexact[];
    scalar pexact[];
    vorticity (u, omega);
    double Luinf = 0, Lpinf = 0, du, dp, maxp = 0, maxu = 0;
    int curlevel=1;
    theory(uexact, pexact, t+dt, RE);
    scalar l[]; foreach(reduction(max:curlevel)) {l[] = level; if (curlevel<level) curlevel=level;}
#ifdef DEBUG_MODE
    vector gradu[], gradv[];
    gradients({u.x, u.y}, {gradu, gradv});
    foreach(){
        mod_du_dx[] = sqrt(sq(gradu.x[]) + sq(gradu.y[]) + sq(gradv.x[]) + sq(gradv.y[]));
    }
#endif
    foreach(reduction(max:Luinf) reduction(max:Lpinf) reduction(max:maxu) reduction(max:maxp)){
        if (fs[] == 0){ //in inner frame
            du = sqrt(sq(u.x[] - uexact.x[]) + sq(u.y[] - uexact.y[]));
            dp = fabs(p[] - pexact[]);
            if (du > Luinf) Luinf = du;
            if (dp > Lpinf) Lpinf = dp;
            if (norm(uexact) > maxu) maxu = norm(uexact);
            if (fabs(pexact[]) > maxp) maxp = fabs(pexact[]);
        }
    }
    int tnc = count_cells();
    int nmax = (int)pow (2, maxlevel*dimension);
    fprintf (ferr, "i= %d t+dt= %g dt= %g Luinf= %g Lpinf= %g Luinf_rel= %g Lpinf_rel= %g eta_s= %g count_cells= %d nmax= %d compress_ratio= %g curlevel= %d\n", i, t+dt, dt, Luinf, Lpinf, Luinf/u0, Lpinf/p0, eta_s, tnc, nmax, (double)tnc/nmax, curlevel);
    double eps_arr[] = {1, 1, 1, 1, 1, 1, 1};
    MinMaxValues((scalar *){p, pexact, u.x, u.y, uexact.x, uexact.y, omega}, eps_arr);
}

//Output
event vtk_file (i++){
//event vtk_file (i+=5){
//event vtk_file (t += 0.01){
    char subname[80]; sprintf(subname, "vortex_BP_periodic_adapt");
    scalar l[]; foreach() l[] = level;
    vector uexact[];
    scalar pexact[];
    double time_of_output = (iter_fp==0) ? 0 : t+dt;
    theory(uexact, pexact, time_of_output, RE);

#ifdef DEBUG_MODE
    double d = 0;
    foreach() {
        d = 0.;
        foreach_dimension(){
            d += uf.x[1] - uf.x[];
        }
        divutmpAfter[] = d/(dt*Delta);
    }
    boundary((scalar*){divutmpAfter});
#endif
    output_vtu_MPI( (scalar *) {fs, omega, p, pf, pexact, l, divutmpAfter, deltap},
            (vector *) {u, uexact, deltag, g, target_Uv, dbp, total_rhs, tmp_err_u},
            (vector *) {uf}, subname,  time_of_output);
//    output_vtu_MPI( (scalar *) {fs, omega, p, pexact, l, pid_num}, (vector *) {u, uexact}, subname);
}
/**
We adapt according to the error on the embedded geometry, velocity*/
//#define ADAPT_SCALARS {omega, mod_du_dx}
//#define ADAPT_EPS_SCALARS (double[]){eps*2*pi, eps*sqrt(2.0)*pi}
event adapt (i++){
//    double eps_arr[] = ADAPT_EPS_SCALARS;//??? sometimes doesn't work
//    MinMaxValues(ADAPT_SCALARS, eps_arr);
//    adapt_wavelet ((scalar *) {fs, omega}, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
    adapt_wavelet ((scalar *) {u.x, u.y}, (double[]){eps, eps}, maxlevel = maxlevel, minlevel = minlevel);
}

//event stop (i = 30){
event stop (t = 0.3){
    event("end_timestep");
}
//event stop (i = 10);
/**
## See also

* [Same example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/cylinder.html)
*/
