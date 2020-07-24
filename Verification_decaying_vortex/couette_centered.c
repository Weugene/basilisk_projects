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

#define DEBUG_MODE
#define FILTERED
//#define JACOBI 1
#define EPS_MAXA 2
//#define RELATIVE_RESIDUAL
#define PRINT_ALL_VALUES
//#define PERIODIC_BC
scalar omega[];
vector target_Uv[];
vector uold[];
vector tmp_err_u[];
vector divtauu[];
scalar my_residual[];
double mydt = 0;
int maxlevel;
int minlevel;
double Ldomain = 1.0, RE = 30., rho1 = 1, G = 1, h = 1, U = 1, mu1;
coord C;
double eps = 1e-2;
double ducor=0;
#ifdef DEBUG_MODE
    scalar divutmp[], divutmpAfter[];
    scalar mod_du_dx[];
    vector conv_term[];
#endif

#include "navier-stokes/centered.h"
#include "../src_local/utils-weugene.h"
#include "../src_local/output_vtu_foreach.h"
#include "fractions.h"
face vector muv[];

//max U = 4.25
#define uxeInCh(x,y,t) ((G/(2.0*mu1))*y*(h-y) + U*y/h)
#define uxe(x,y,t) uxeInCh(x,y,t)
#define uye(x,y,t) 0
#define pe(x,y,t)  G*(0.5*Ldomain - x)
#define transientTime 1e+5
#define funtran(t) min((t)/transientTime, 1)
//#define uxeBC(x,y,t) (uxe(x,y,t))
#define uxeBC(x,y,t) ((y)*(1-funtran(t)) + funtran(t)*uxe(x,y,t))
#define uyeBC(x,y,t) (uye(x,y,t))
#define outBC(var,Delta_ghost,dt_) (var[] + ((var[] - var[-1])/Delta)*(Delta_ghost-u.x[]*dt_)/(1.0 + dt_*(u.x[] - u.x[-1])/Delta))// set on ghost cells
//#define outxBC(x_,y_,t_,dt_) (Delta*(-(u.x[] - uold.x[])/(dt_*C.x)) + u.x[])
//#define outyBC(x_,y_,t_,dt_) (Delta*(-(u.y[] - uold.y[])/(dt_*C.x)) + u.y[])
//#define outxBC(x_,y_,t_,dt_) (Delta*(-(u.x[] - uold.x[])/(dt_*u.x[])) + u.x[])
//#define outyBC(x_,y_,t_,dt_) (Delta*(-(u.y[] - uold.y[])/(dt_*u.x[])) + u.y[])
//#define outxBC(x_,y_,t_,dt_) neumann(-(u.x[] - uold.x[])/(dt_*C.x))
//#define outyBC(x_,y_,t_,dt_) neumann(-(u.y[] - uold.y[])/(dt_*C.x))

#define u0 1
#define p0 1

#define dpdx(x_,y_,t_) (-G)
#define dpdy(x_,y_,t_) 0

#define dudxe 0
#define dudye 0

void theory(vector u, scalar p, double t, double RE){
        foreach(){
            u.x[] = uxe(x,y,t);
            u.y[] = uye(x,y,t);
            p[]   = pe(x,y,t);
        }
        boundary({u, p});
        fprintf(ferr, "theory  t= %g\n", t);
}

void calculateC(){
    double Cx=0, Cy=0;
    foreach_boundary(left, reduction(+:Cx) reduction(+:Cy)){
        Cx += uf.x[]*Delta;
        Cy += uf.y[-1,0]*Delta;
    }
    C.x = Cx/Ldomain;
    C.y = Cy/Ldomain;
    fprintf(ferr, "mean mass velocity at left: C=%g %g\n", C.x, C.y);

    Cx=0, Cy=0;
    foreach_boundary(right, reduction(+:Cx) reduction(+:Cy)){
        Cx += uf.x[ghost]*Delta;
        Cy += uf.y[ghost]*Delta;
    }
    Cx /= Ldomain;
    Cy /= Ldomain;
    ducor = C.x - Cx;
    fprintf(ferr, "mean mass velocity at right: C=%g %g\n", Cx, Cy);
    fprintf(ferr, "difference: deltaC=%g %g\n", C.x - Cx, C.y - Cy);
}
void deepCopy(vector u, vector uold){
    foreach() foreach_dimension() uold.x[] = u.x[];
    for (int b = 0; b < nboundary; b++) foreach_boundary(b) foreach_dimension() {
        uold.x[ghost] = u.x[ghost];
    }
}

//in 2D each cells must have 2 BC, in 3D - 3BC

// u^{n+1} = u^n + F(x,y,t,C,u^n)*dt
u.n[left]   = dirichlet(uxeBC(x,y,t)); // available x,y,z at the computational domain
u.n[right]  = dirichlet(uxeBC(x,y,t));
u.n[bottom] = dirichlet(0);
u.n[top]    = dirichlet(0);

u.t[left]   = dirichlet(0);
u.t[right]  = dirichlet(0);
u.t[bottom] = dirichlet(0);
u.t[top]    = dirichlet(1);

//uf.n[left]   = uxeBC(x,y,t);
//uf.n[bottom] = 0;
//uf.n[top]    = 0;
//
//uf.t[left]   = 0;
//uf.t[bottom] = 0;
//uf.t[top]    = 1;

//p[left]   = neumann(rho1*(uf.x[] -      uxeBC(x,y,t+dt))/dt);
//p[right]  = neumann(rho1*(uf.x[] -      outBC(uf.x,0.5*Delta,t+dt) - ducor)/dt);
//p[bottom] = neumann(rho1*(uf.y[] -      0)/dt);
//p[top]    = neumann(rho1*(uf.y[ghost] - 0)/dt);
//
//pf[left]   = neumann(rho1*(uf.x[] -      uxeBC(x,y,t+0.5*dt))/(0.5*dt));
//pf[right]  = neumann(rho1*(uf.x[ghost] - outBC(uf.x,0.5*Delta,t+0.5*dt) - ducor)/(0.5*dt));
//pf[bottom] = neumann(rho1*(uf.y[] -      0)/(0.5*dt));
//pf[top]    = neumann(rho1*(uf.y[ghost] - 0)/(0.5*dt));

int main(int argc, char * argv[]){
        minlevel = 4;
        maxlevel = 7;
//        eta_s=1e-5;
//    cells_per_zone = 2;
        eps = 1e-3;
#ifdef PERIODIC_BC
        periodic(right);
#endif
        size (Ldomain);
        origin (-0.5*Ldomain, 0);
        DT = 1e-3;
        CFL = 0.1;
        TOLERANCE = 6e-10;
        RELATIVE_RES_TOLERANCE = 0.05;
        NITERMAX = 30;
        reference_pressure = true;
        mu1 = 1./RE;
        mu = muv;
        C.x = (h/Ldomain)*(G*sq(h) + 6*U*mu1)/(12.*mu1);
        C.y = 0;
        N = 1<<maxlevel;
        fprintf(ferr, "maxlevel= %d eps= %g\n", maxlevel, eps);
        fprintf(ferr, "TOL=%g NITERMAX=%d Re=%g Rel_res_tol=%g C=%g %g\n", TOLERANCE, NITERMAX, RE, RELATIVE_RES_TOLERANCE, C.x, C.y);
        fprintf(ferr, "CFL=%g eps=%g\n", CFL, eps);
        run();
}


event init (t = 0){
        vector uexact[];
        scalar pexact[];
//        const face vector mdpdx[] = {1, 0};
//        a = mdpdx;
        if (!restore (file = "restart")) {
            int it = 0;
//        do {
            event ("properties");
            theory(uexact, pexact, 0, RE);
            foreach() {
                u.x[] = y;
                u.y[] = 0;
            }
            boundary({u});
//        } while ( ++it <= 10 && adapt_wavelet((scalar *){fs, u}, (double []){eps, eps, eps}, maxlevel=maxlevel, minlevel=minlevel).nf != 0);
//            foreach() {
//                g.x[] = 1;
//                p[] = 0;// my initial guess
//                pf[] = 0;// my initial guess
//            }
            deepCopy(u, uold);
        }
//        event("vtk_file");
}

/**
We set a constant viscosity corresponding to a Reynolds number of 40, 100,
based on the cylinder diameter (1) and the inflow velocity (1). */

event properties (i++){
    foreach_face() muv.x[] = mu1;
    boundary((scalar *){muv});
}

event set_dtmax (i++){
    NITERMIN = 1;
    NITERMAX = 100;
    DT *= 1.05;
    DT = min(DT, 0.001);
    fprintf(ferr, "set_dtmax: tnext= %g Dt= %g", tnext, DT);
}



event advection_term(i++){
    calculateC ();
    deepCopy(u, uold);
    u.n[left]  = dirichlet(uxeBC(x,y,t+dt));
    u.n[right] = outBC(u.x,Delta,dt);
    u.t[right] = outBC(u.y,Delta,dt);

//    uf.n[left]   = uxeBC(x,y,t+0.5*dt);
//    uf.n[right]  = outBC(uf.x,0.5*Delta,0.5*dt);
//    uf.n[bottom] = 0;
//    uf.n[top]    = 0;
//
//    uf.t[left]   = 0;
//    uf.t[right]  = outBC(uf.y,0.5*Delta,0.5*dt);
//    uf.t[bottom] = 0;
//    uf.t[top]    = 1;
}


event projection(i++){
    calculateC ();
//    uf.n[left]   = uxeBC(x,y,t+dt);
//    uf.n[right]  = outBC(uf.x,0.5*Delta,dt);
//    uf.n[bottom] = 0;
//    uf.n[top]    = 0;
//
//    uf.t[left]   = 0;
//    uf.t[right]  = outBC(uf.y,0.5*Delta,dt);
//    uf.t[bottom] = 0;
//    uf.t[top]    = 1;
}
//event end_timestep (i += 1){
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
        du = sqrt(sq(u.x[] - uexact.x[]) + sq(u.y[] - uexact.y[]));
        dp = fabs(p[] - pexact[]);
        if (du > Luinf) Luinf = du;
        if (dp > Lpinf) Lpinf = dp;
        if (norm(uexact) > maxu) maxu = norm(uexact);
        if (fabs(pexact[]) > maxp) maxp = fabs(pexact[]);
    }

    int tnc = count_cells();
    int nmax = (int)pow (2, maxlevel*dimension);
    fprintf (ferr, "i= %d t+dt= %g dt= %g Luinf= %g Lpinf= %g Luinf_rel= %g Lpinf_rel= %g count_cells= %d nmax= %d compress_ratio= %g curlevel= %d\n", i, t+dt, dt, Luinf, Lpinf, Luinf/u0, Lpinf/p0, tnc, nmax, (double)tnc/nmax, curlevel);
//    double eps_arr[] = {1, 1, 1, 1, 1, 1, 1};
//    MinMaxValues((scalar *){p, pexact, u.x, u.y, uexact.x, uexact.y, omega}, eps_arr);
}

//Output
//event vtk_file (i++){
//event vtk_file (i+=10){
event vtk_file (t += 1){
    char subname[80]; sprintf(subname, "couette_centered2");
    scalar l[]; foreach() {l[] = level;}
    vector uexact[];
    scalar pexact[];
    theory(uexact, pexact, t+mydt, RE);

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
    output_vtu_MPI( (scalar *) {omega, p, pf, pexact, l, divutmpAfter, divutmp, my_residual},
    (vector *) {u, uold, uexact, g, divtauu},
    (vector *) {uf}, subname, t + mydt );
    fprintf(ferr, "dt=%g mydt=%g\n", dt, mydt);
}
/**
We adapt according to the error on the embedded geometry, velocity*/
//#define ADAPT_SCALARS {omega, mod_du_dx}
//#define ADAPT_EPS_SCALARS (double[]){eps*2*pi, eps*sqrt(2.0)*pi}
//event adapt (i++){
////    double eps_arr[] = ADAPT_EPS_SCALARS;//??? sometimes doesn't work
////    MinMaxValues(ADAPT_SCALARS, eps_arr);
////    adapt_wavelet ((scalar *) {fs, omega}, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
//    adapt_wavelet ((scalar *) {u.x, u.y}, (double[]){eps, eps}, maxlevel = maxlevel, minlevel = minlevel);
//}

//event stop (i = 10){
event stop (t = 50){
    event("end_timestep");
};

/**
## See also

* [Same example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/cylinder.html)
*/