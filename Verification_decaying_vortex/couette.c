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
//#define DEBUG_BRINKMAN_PENALIZATION 1
//#define DEBUG_OUTPUT_VTU_MPI
#define DEBUG_MINMAXVALUES
#define DEBUG_MODE
#define FILTERED
#define PRINT_ALL_VALUES
//#define JACOBI 1
#define EPS_MAXA 2
#define RELATIVE_RESIDUAL
#define MODIFIED_CHORIN 0
#define PERIODIC_BC
scalar omega[], fs[];
vector target_Uv[];
vector deltag[];
scalar deltap[];
scalar phi[];
vector uold[];
vector tmp_err_u[];
vector divtauu[];
face vector ufold[];
double mydt = 0;
int maxlevel;
int minlevel;
double Ldomain = 1.0, RE = 30., rho1 = 1, G = 1, h = 1, U = 1, mu1;
coord C;
double eps = 1e-2;
#ifdef DEBUG_MODE
    scalar divutmp[], divutmpAfter[];
    scalar mod_du_dx[];
    vector conv_term[];
#endif

#include "../src_local/centered-weugene.h"
#include "../src_local/output_vtu_foreach.h"
#include "fractions.h"
face vector muv[];

void frame_BP(scalar fs, double t){
    vertex scalar phi[];
    face vector fs_face[];
    foreach_vertex() {
//        phi[] = 0;
//        phi[] = -( pow(x,10) + pow(y,10) - 0.117 );
//        phi[] = ( fabs(x) <= 1 && fabs(y) <= 1 ) ? 0 : 1;
        phi[] = ( y*(y-h) + Ldomain/pow(2, maxlevel+1) );
    }
    boundary ({phi});
    fractions (phi, fs, fs_face);
}

//max U = 4.25
#define uxeInCh(x,y,t) ((G/(2.0*mu1))*y*(h-y) + U*y/h)
#define uxe(x,y,t) uxeInCh(x,y,t)
#define uye(x,y,t) 0
#define pe(x,y,t)  G*(0.5*Ldomain - x)
#define transientTime 5
#define uxeBC(x,y,t) dirichlet(sq(y)*max((transientTime-t)/transientTime, 0) + min(t/transientTime, 1)*uxe(x,y,t))
#define uyeBC(x,y,t) dirichlet(uye(x,y,t))
#define outxBC(x_,y_,t_,dt_) neumann(-(u.x[] - uold.x[])/(dt_*u.x[]))
#define outyBC(x_,y_,t_,dt_) neumann(-(u.y[] - uold.y[])/(dt_*u.x[]))

//#define outxBC(x_,y_,t_,dt_) neumann((uold.x[] + uold.x[ghost])/(2.0*dt_*uold.x[] - Delta))*(2.0*dt_*uold.x[]/Delta - 1)/(2.0*dt_*uold.x[]/Delta + 1)
//#define outyBC(x_,y_,t_,dt_) neumann((uold.y[] + uold.y[ghost])/(2.0*dt_*uold.x[] - Delta))*(2.0*dt_*uold.x[]/Delta - 1)/(2.0*dt_*uold.x[]/Delta + 1)
//#define outxBC(x_,y_,t_,dt_) (uold.x[] + uold.x[ghost] + (2.0*(dt_*uold.x[]/Delta) - 1)*u.x[])/(1.0 + 2.0*dt_*uold.x[]/Delta)
//#define outyBC(x_,y_,t_,dt_) (uold.y[] + uold.y[ghost] + (2.0*(dt_*uold.x[]/Delta) - 1)*u.y[])/(1.0 + 2.0*dt_*uold.x[]/Delta)
#define outxBCface(x_,y_,t_,dt_) neumann(0)
#define outyBCface(x_,y_,t_,dt_) neumann(0)
//#define outxBCface(x_,y_,t_,dt_) neumann(-(u.x[] - uold.x[])/(dt_*C.x))
//#define outyBCface(x_,y_,t_,dt_) neumann(-(u.y[] - uold.y[])/(dt_*C.x))
//#define outxBCface(x_,y_,t_,dt_) ((ufold.x[ghost] + (uold.x[]*dt_/Delta)*uf.x[])/(1 + (uold.x[]*dt_/Delta)))
//#define outyBCface(x_,y_,t_,dt_) ((ufold.y[ghost] + (uold.x[]*dt_/Delta)*uf.y[])/(1 + (uold.x[]*dt_/Delta)))
//#define uxeBC(x,y,t) neumann(0)
//#define uyeBC(x,y,t) neumann(0)
//#define uxeBC(x,y,t) (fs[]>0) ? neumann(0) : dirichlet(uxe(x,y,t))
//#define uyeBC(x,y,t) (fs[]>0) ? neumann(0) : dirichlet(uye(x,y,t))
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
//in 2D each cells must have 2 BC, in 3D - 3BC

// u^{n+1} = u^n + F(x,y,t,C,u^n)*dt
//u.n[left]   = dirichlet(uxe(x,y,t)); // available x,y,z at the computational domain
    u.n[left]   = uxeBC(x,y,t); // available x,y,z at the computational domain
    u.n[right]  = uxeBC(x,y,t);
    u.n[bottom] = dirichlet(0);
    u.n[top]    = dirichlet(0);

//u.t[left]   = dirichlet(0);
    u.t[left]   = dirichlet(0);
    u.t[right]  = dirichlet(0);
    u.t[bottom] = dirichlet(0);
    u.t[top]    = dirichlet(1);

//uold.n[left]   = dirichlet(uxe(x,y,t)); // available x,y,z at the computational domain
    uold.n[left]   = uxeBC(x,y,t); // available x,y,z at the computational domain
    uold.n[right]  = uxeBC(x,y,t);
    uold.n[bottom] = dirichlet(0);
    uold.n[top]    = dirichlet(0);

//uold.t[left]   = dirichlet(0);
    uold.t[left]   = dirichlet(0);
    uold.t[right]  = dirichlet(0);
    uold.t[bottom] = dirichlet(0);
    uold.t[top]    = dirichlet(1);

    uf.n[left]   = uxe(x,y,t);
    uf.n[right]  = uxe(x,y,t);
    uf.n[bottom] = 0;
    uf.n[top]    = 0;

    uf.t[left]   = 0;
    uf.t[right]  = 0;
    uf.t[bottom] = 0;
    uf.t[top]    = 1;

//    ufold.n[left]   = uxe(x,y,t);
//    ufold.n[right]  = uxe(x,y,t);
//    ufold.n[bottom] = 0;
//    ufold.n[top]    = 0;
//
//    ufold.t[left]   = 0;
//    ufold.t[right]  = 0;
//    ufold.t[bottom] = 0;
//    ufold.t[top]    = 1;
//u.t[right] = neumann(0);
//u.n[right] = neumann(0);
//u.n[right] = neumann(-u.x[ghost]/(C*dt))/(1-Delta/(C*dt));
//p[right]   = neumann(-G);

//uf.n[right]  = neumann(0);
//uf.t[right]  = neumann(0);
//uf.n[right]  = dirichlet(uxe(x,y,t+dt));
//uf.t[right]  = dirichlet(uxe(x,y,t+dt));

//p[left]   = neumann(G);
//p[right]  = dirichlet(pe(x,y,t));
//p[bottom] = dirichlet(pe(x,y,t));
//p[top]    = dirichlet(pe(x,y,t));

//p[left]   = neumann(G);
//p[right]  = neumann(G);
//p[bottom] = neumann(0);
//p[top]    = neumann(0);

    p[left]   = neumann(rho1*(uf.x[] -      uxe(x,y,t+dt))/dt);
    p[right]  = neumann(0);
    p[bottom] = neumann(rho1*(uf.y[] -      uye(x,y,t+dt))/dt);
    p[top]    = neumann(rho1*(uf.y[ghost] - uye(x,y,t+dt))/dt);

    pf[left]   = neumann(rho1*(uf.x[]      - uxe(x,y,t + 0.5*dt))/(0.5*dt));
    pf[right]  = neumann(0);
    pf[bottom] = neumann(rho1*(uf.y[]      - uye(x,y,t + 0.5*dt))/(0.5*dt));
    pf[top]    = neumann(rho1*(uf.y[ghost] - uye(x,y,t + 0.5*dt))/(0.5*dt));

    deltap[left]   = neumann(rho1*(uf.x[]      - uxe(x,y,t+dt))/dt);
    deltap[right]  = neumann(0);
    deltap[bottom] = neumann(rho1*(uf.y[]      - uye(x,y,t+dt))/dt);
    deltap[top]    = neumann(rho1*(uf.y[ghost] - uye(x,y,t+dt))/dt);


//deltap[left]   = neumann(0);
//deltap[right]  = neumann(0);
//deltap[bottom] = neumann(0);
//deltap[top]    = neumann(0);

    phi[left]   = neumann(uf.x[]      - uxe(x,y,t));
    phi[right]  = neumann(uf.x[ghost] - uxe(x,y,t));
    phi[bottom] = neumann(uf.y[]      - uye(x,y,t));
    phi[top]    = neumann(uf.y[ghost] - uye(x,y,t));

    int main(int argc, char * argv[]) {
        minlevel = 4;
        maxlevel = 7;
//        eta_s=1e-5;
//    cells_per_zone = 2;
        eps = 1e-3;
        if (argc > 1) {
            maxlevel = atoi(argv[1]);
        }
//        if (argc > 2) {
//            eta_s = atof(argv[2]);
//        }
        if (argc > 3) {
            eps = atof(argv[3]);
        }
#ifdef PERIODIC_BC
        periodic(right);
#endif
        size (Ldomain);
        origin (-0.5*Ldomain, 0);
        DT = 1e-2;
        CFL = 0.1;
        TOLERANCE = 1e-8;
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

    void update_targetU(vector target_Uv, double t){
        foreach() {
            target_Uv.x[] = (y>=0.5)? 1 : 0;
            target_Uv.y[] = 0;
        }
        boundary((scalar *){target_Uv});
    }

    event init (t = 0)
    {
        vector uexact[];
        scalar pexact[];
        if (!restore (file = "restart")) {
            int it = 0;
//        do {
            event ("properties");
            theory(uexact, pexact, 0, RE);
            foreach() {
                u.x[] = uexact.x[];
                u.y[] = 0;
            }
//            foreach() u.x[] = 0;
            boundary({u});
//        } while ( ++it <= 10 && adapt_wavelet((scalar *){fs, u}, (double []){eps, eps, eps}, maxlevel=maxlevel, minlevel=minlevel).nf != 0);
            foreach() {
                g.x[] = 1;
                p[] = 0;// my initial guess
                pf[] = 0;// my initial guess
                foreach_dimension() uold.x[] = u.x[];
            }
            boundary({p, pf, uold, g});
        }
        event("vtk_file");
    }

/**
We set a constant viscosity corresponding to a Reynolds number of 40, 100,
based on the cylinder diameter (1) and the inflow velocity (1). */

    event properties (i++){

    foreach_face() muv.x[] = mu1;
    boundary((scalar *){muv});
    }

    event set_dtmax (i++){
    NITERMIN = 5;
    NITERMAX = 100;
    DT *= 1.05;
    DT = min(DT, 0.01);
    fprintf(ferr, "set_dtmax: tnext= %g Dt= %g", tnext, DT);
    }

    event advection_term(i = 0){
//    u.n[right] = neumann(0);
//    u.t[right] = neumann(0);
//    uold.n[right] = neumann(0);
//    uold.t[right] = neumann(0);
        u.n[right] = outxBC(x,y,t,dt);
        u.t[right] = outyBC(x,y,t,dt);
//        uold.n[right] = outxBC(x,y,t,dt);
//        uold.t[right] = outyBC(x,y,t,dt);

//        uf.n[right]  = outxBCface(x,y,t,dt);
//        uf.t[right]  = outyBCface(x,y,t,dt);
//        ufold.n[right]  = outxBCface(x,y,t,dt);
//        ufold.t[right]  = outyBCface(x,y,t,dt);
    }

    event advection_term(i++){
//    uf.n[right]  = outxBC(x,y,t,0.5*dt);
//    uf.t[right]  = outyBC(x,y,t,0.5*dt);
//    uf.n[right]  = neumann(0);
//    uf.t[right]  = neumann(0);
    //deep copy of uold with ghost cells
    uf.n[right]  = outxBCface(x,y,t,0.5*dt);
    uf.t[right]  = outyBCface(x,y,t,0.5*dt);
//    ufold.n[right]  = outxBCface(x,y,t,0.5*dt);
//    ufold.t[right]  = outyBCface(x,y,t,0.5*dt);
    foreach() foreach_dimension() uold.x[] = u.x[];
    foreach_face() {
        ufold.x[] = uf.x[];
    }
//    boundary({ufold});
    for (int b = 0; b < nboundary; b++) foreach_boundary(b) foreach_dimension() {
        uold.x[ghost] = u.x[ghost];
    }
    }

event projection(i++){
//    foreach_face() {
//        ufold.x[] = uf.x[];
//    }
    uf.n[right]    = outxBCface(x,y,t,dt);
    uf.t[right]    = outyBCface(x,y,t,dt);
//    ufold.n[right] = outxBCface(x,y,t,dt);
//    ufold.t[right] = outyBCface(x,y,t,dt);
//    uf.n[right]  = outxBC(x,y,t,dt);
//    uf.t[right]  = outyBC(x,y,t,dt);
//    uf.n[right]  = neumann(0);
//    uf.t[right]  = neumann(0);
}

//event end_timestep(i++){
//    double xref = -1.5, yref = -1.5, zref = 0; //reference location
//    double pref = 0;                                            //reference pressure
//    double pcor = interpolate (p, xref, yref, zref) - pref;
//    fprintf(ferr, "pressure correction=%g at x=%g y=%g z=%g\n", pcor, xref, yref, zref);
//#if _MPI
//    MPI_Allreduce (MPI_IN_PLACE, &pcor, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
//#endif
//    foreach() p[] -= pcor;
//    boundary ({p});
//}

event end_timestep(i++){
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
    Cx += uf.x[1,0]*Delta;
    Cy += uf.y[1,0]*Delta;
    }
    fprintf(ferr, "mean mass velocity at right: C=%g %g\n", Cx/Ldomain, Cy/Ldomain);
    }
event end_timestep (i += 1){
//event end_timestep (i += 100){
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
    fprintf (ferr, "i= %d t+dt= %g dt= %g Luinf= %g Lpinf= %g Luinf_rel= %g Lpinf_rel= %g count_cells= %d nmax= %d compress_ratio= %g curlevel= %d\n", i, t+dt, dt, Luinf, Lpinf, Luinf/u0, Lpinf/p0, tnc, nmax, (double)tnc/nmax, curlevel);
    double eps_arr[] = {1, 1, 1, 1, 1, 1, 1};
    MinMaxValues((scalar *){p, pexact, u.x, u.y, uexact.x, uexact.y, omega}, eps_arr);
    }

//Output
event vtk_file (i++){
//event vtk_file (i+=10){
//event vtk_file (t += 1){
    char subname[80]; sprintf(subname, "couette_not_periodic_adapt_jac2");
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
    output_vtu_MPI( (scalar *) {phi, omega, p, pf, pexact, l, divutmpAfter, deltap, divutmp},
    (vector *) {u, uold, uexact, deltag, g, divtauu},
    (vector *) {uf, ufold}, subname, t + mydt );
//    (vector *) {u, uexact, uf_low, uf_up, deltag, g, target_Uv, u_low, u_up, p_low, p_up, dbp, total_rhs, tmp_err_u, divtauu}, subname, t + mydt );
//    output_vtu_MPI( (scalar *) {fs, omega, p, pexact, l, pid_num}, (vector *) {u, uexact}, subname, periodic_BC);
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




//    /**
//# Decaying vortex problem
//To verify the spatial accuracy of the present method, the decaying vortex problem is chosen because it is an unsteady problem with an analytical solution:
// $$u(x, y, t) = − \cos 􏱲x \sin 􏱲y \exp{−2􏱲\pi^2 t/Re}$$
// $$v(x, y, t) =   \sin 􏱲x \cos 􏱲y \exp{−2􏱲\pi^2 t/Re}$$
// $$p(x, y, t) =  −\frac14(\cos 2􏱲x + \cos 2􏱲y)\exp{−4\pi^2 t/Re}$$
//
//The computational domain is −1.5 < x, y 􏰒< 1.5 and the IB is located at x = ± 1 and y = ± 1.
//The Reynolds number based on the maximum velocity and vortex size is set to 30, and the initial and boundary conditions are given by the analytical solution above.
//Simulations are performed till $t=0.3$.
//We use the centered Navier-Stokes solver, with embedded boundaries and
//advect the passive tracer *f*. */
//
//#define BRINKMAN_PENALIZATION 1
//#define DEBUG_BRINKMAN_PENALIZATION 1
////#define DEBUG_OUTPUT_VTU_MPI
//#define DEBUG_MINMAXVALUES
//#define DEBUG_MODE
//#define FILTERED
////#define JACOBI 1
//#define EPS_MAXA 2
//#define RELATIVE_RESIDUAL
//#define MODIFIED_CHORIN 0
////#define PERIODIC_BC
//    scalar omega[], fs[];
//    vector target_Uv[];
//    vector deltag[];
//    scalar deltap[];
//    scalar phi[];
//    vector uold[];
//    vector tmp_err_u[];
//    vector divtauu[];
//
//    double mydt = 0;
//    int maxlevel;
//    int minlevel;
//    double Ldomain = 1.0, RE = 30., rho1 = 1, G = 1, h = 1, U = 1, mu1;
//    coord C;
//    double eps = 1e-2;
//#ifdef DEBUG_MODE
//    scalar divutmp[], divutmpAfter[];
//    scalar mod_du_dx[];
//    vector conv_term[];
//#endif
//
//#include "../src_local/centered-weugene.h"
////#include "view.h"
//#include "../src_local/output_vtu_foreach.h"
//#include "fractions.h"
//    face vector muv[];
//
//    void frame_BP(scalar fs, double t){
//        vertex scalar phi[];
//        face vector fs_face[];
//        foreach_vertex() {
////        phi[] = 0;
////        phi[] = -( pow(x,10) + pow(y,10) - 0.117 );
////        phi[] = ( fabs(x) <= 1 && fabs(y) <= 1 ) ? 0 : 1;
//            phi[] = ( y*(y-h) + Ldomain/pow(2, maxlevel+1) );
//        }
//        boundary ({phi});
//        fractions (phi, fs, fs_face);
//    }
//
////max U = 4.25
//#define uxeInCh(x,y,t) ((G/(2.0*mu1))*y*(h-y) + U*y/h)
//#define uxe(x,y,t) ((1.0 - fs[])*uxeInCh(x,y,t) + fs[]*target_Uv.x[])
//#define uye(x,y,t) 0
//#define pe(x,y,t)  G*(0.5*Ldomain - x)
//#define uxeBC(x,y,t) dirichlet(uxe(x,y,t))
//#define uyeBC(x,y,t) dirichlet(uye(x,y,t))
//#define outxBC(x_,y_,t_,dt_) neumann((uold.x[] + uold.x[ghost])/(2.0*dt_*C.x - Delta))*(2.0*dt_*C.x/Delta - 1)/(2.0*dt_*C.x/Delta + 1)
//#define outyBC(x_,y_,t_,dt_) neumann((uold.y[] + uold.y[ghost])/(2.0*dt_*C.y - Delta))*(2.0*dt_*C.y/Delta - 1)/(2.0*dt_*C.y/Delta + 1)
////#define outxBC(x_,y_,t_,dt_) (uold.x[] + uold.x[ghost] + (2.0*(dt*C/Delta) - 1)*u.x[])/(1.0 + 2.0*dt*C/Delta)
////#define outyBC(x_,y_,t_,dt_) (uold.y[] + uold.y[ghost] + (2.0*(dt*C/Delta) - 1)*u.y[])/(1.0 + 2.0*dt*C/Delta)
//#define outxBCface(x_,y_,t_,dt_)
//#define outyBCface(x_,y_,t_,dt_)
////#define uxeBC(x,y,t) neumann(0)
////#define uyeBC(x,y,t) neumann(0)
////#define uxeBC(x,y,t) (fs[]>0) ? neumann(0) : dirichlet(uxe(x,y,t))
////#define uyeBC(x,y,t) (fs[]>0) ? neumann(0) : dirichlet(uye(x,y,t))
//#define u0 1
//#define p0 1
//
//#define dpdx(x_,y_,t_) (-G)
//#define dpdy(x_,y_,t_) 0
//
//#define dudxe 0
//#define dudye 0
//    void theory(vector u, scalar p, double t, double RE){
//        foreach(){
//            u.x[] = uxe(x,y,t);
//            u.y[] = uye(x,y,t);
//            p[]   = pe(x,y,t);
//        }
//        boundary({u, p});
//        fprintf(ferr, "theory  t= %g\n", t);
//    }
////in 2D each cells must have 2 BC, in 3D - 3BC
//
//// u^{n+1} = u^n + F(x,y,t,C,u^n)*dt
////u.n[left]   = dirichlet(uxe(x,y,t)); // available x,y,z at the computational domain
//    u.n[left]   = uxeBC(x,y,t); // available x,y,z at the computational domain
//    u.n[right]  = uxeBC(x,y,t);
//    u.n[bottom] = uyeBC(x,y,t);
//    u.n[top]    = uyeBC(x,y,t);
//
////u.t[left]   = dirichlet(0);
//    u.t[left]   = uyeBC(x,y,t);
//    u.t[right]  = uyeBC(x,y,t);
//    u.t[bottom] = uxeBC(x,y,t);
//    u.t[top]    = uxeBC(x,y,t);
//
////uold.n[left]   = dirichlet(uxe(x,y,t)); // available x,y,z at the computational domain
//    uold.n[left]   = uxeBC(x,y,t); // available x,y,z at the computational domain
//    uold.n[right]  = uxeBC(x,y,t);
//    uold.n[bottom] = uyeBC(x,y,t);
//    uold.n[top]    = uyeBC(x,y,t);
//
////uold.t[left]   = dirichlet(0);
//    uold.t[left]   = uyeBC(x,y,t);
//    uold.t[right]  = uyeBC(x,y,t);
//    uold.t[bottom] = uxeBC(x,y,t);
//    uold.t[top]    = uxeBC(x,y,t);
//
//    uf.n[left]   = uxe(x,y,t);
//    uf.n[right]  = uxe(x,y,t);
//    uf.n[bottom] = uye(x,y,t);
//    uf.n[top]    = uye(x,y,t);
//
//    uf.t[left]   = uye(x,y,t);
//    uf.t[right]  = uye(x,y,t);
//    uf.t[bottom] = uxe(x,y,t);
//    uf.t[top]    = uxe(x,y,t);
//
////u.t[right] = neumann(0);
////u.n[right] = neumann(0);
////u.n[right] = neumann(-u.x[ghost]/(C*dt))/(1-Delta/(C*dt));
////p[right]   = neumann(-G);
//
////uf.n[right]  = neumann(0);
////uf.t[right]  = neumann(0);
////uf.n[right]  = dirichlet(uxe(x,y,t+dt));
////uf.t[right]  = dirichlet(uxe(x,y,t+dt));
//
////p[left]   = neumann(G);
////p[right]  = dirichlet(pe(x,y,t));
////p[bottom] = dirichlet(pe(x,y,t));
////p[top]    = dirichlet(pe(x,y,t));
//
////p[left]   = neumann(G);
////p[right]  = neumann(G);
////p[bottom] = neumann(0);
////p[top]    = neumann(0);
//
//    p[left]   = neumann(rho1*(uf.x[] -      uxe(x,y,t+dt))/dt);
//    p[right]  = neumann(rho1*(uf.x[ghost] - uxe(x,y,t+dt))/dt);
//    p[bottom] = neumann(rho1*(uf.y[] -      uye(x,y,t+dt))/dt);
//    p[top]    = neumann(rho1*(uf.y[ghost] - uye(x,y,t+dt))/dt);
//
//    pf[left]   = neumann(rho1*(uf.x[]      - uxe(x,y,t + 0.5*dt))/(0.5*dt));
//    pf[right]  = neumann(rho1*(uf.x[ghost] - uxe(x,y,t + 0.5*dt))/(0.5*dt));
//    pf[bottom] = neumann(rho1*(uf.y[]      - uye(x,y,t + 0.5*dt))/(0.5*dt));
//    pf[top]    = neumann(rho1*(uf.y[ghost] - uye(x,y,t + 0.5*dt))/(0.5*dt));
//
//    deltap[left]   = neumann(rho1*(uf.x[]      - uxe(x,y,t+dt))/dt);
//    deltap[right]  = neumann(rho1*(uf.x[ghost] - uxe(x,y,t+dt))/dt);
//    deltap[bottom] = neumann(rho1*(uf.y[]      - uye(x,y,t+dt))/dt);
//    deltap[top]    = neumann(rho1*(uf.y[ghost] - uye(x,y,t+dt))/dt);
//
//
////deltap[left]   = neumann(0);
////deltap[right]  = neumann(0);
////deltap[bottom] = neumann(0);
////deltap[top]    = neumann(0);
//
//    phi[left]   = neumann(uf.x[]      - uxe(x,y,t));
//    phi[right]  = neumann(uf.x[ghost] - uxe(x,y,t));
//    phi[bottom] = neumann(uf.y[]      - uye(x,y,t));
//    phi[top]    = neumann(uf.y[ghost] - uye(x,y,t));
//
//    int main(int argc, char * argv[]) {
//        minlevel = 4;
//        maxlevel = 7;
//        eta_s=1e-5;
////    cells_per_zone = 2;
//        eps = 1e-3;
//        if (argc > 1) {
//            maxlevel = atoi(argv[1]);
//        }
//        if (argc > 2) {
//            eta_s = atof(argv[2]);
//        }
//        if (argc > 3) {
//            eps = atof(argv[3]);
//        }
//#ifdef PERIODIC_BC
//        periodic(right);
//#endif
//        size (Ldomain);
//        origin (-0.5*Ldomain, 0);
//        DT = 1e-2;
//        CFL = 0.1;
//        TOLERANCE = 1e-8;
//        RELATIVE_RES_TOLERANCE = 0.05;
//        NITERMAX = 30;
//        reference_pressure = true;
//        mu1 = 1./RE;
//        mu = muv;
//        C.x = 0.25*(G*sq(h) + 6*U*mu1)/(12.*mu1) + 0.25*U;
//        C.y = 0;
//        target_U = target_Uv;
//        N = 1<<maxlevel;
//        fprintf(ferr, "maxlevel= %d eta= %g eps= %g\n", maxlevel, eta_s, eps);
//        fprintf(ferr, "TOL=%g NITERMAX=%d Re=%g Rel_res_tol=%g C=%g %g\n", TOLERANCE, NITERMAX, RE, RELATIVE_RES_TOLERANCE, C.x, C.y);
//        fprintf(ferr, "CFL=%g eps=%g\n", CFL, eps);
//        run();
//    }
//
//    void update_targetU(vector target_Uv, double t){
//        foreach() {
//            target_Uv.x[] = (y>=0.5)? 1 : 0;
//            target_Uv.y[] = 0;
//        }
//        boundary((scalar *){target_Uv});
//    }
//
//    event init (t = 0)
//    {
//        vector uexact[];
//        scalar pexact[];
//        if (!restore (file = "restart")) {
//            int it = 0;
////        do {
//            event ("properties");
//            theory(uexact, pexact, 0, RE);
////            foreach() foreach_dimension() u.x[] = uexact.x[];
//            foreach() u.x[] = 0;
//            boundary({u});
////        } while ( ++it <= 10 && adapt_wavelet((scalar *){fs, u}, (double []){eps, eps, eps}, maxlevel=maxlevel, minlevel=minlevel).nf != 0);
//            foreach() {
//                g.x[] = 1;
//                p[] = 0;//pexact[];// my initial guess
//                pf[] = 0;//pexact[];// my initial guess
//                foreach_dimension() uold.x[] = u.x[];
//            }
//            boundary({p, pf, uold, g});
//        }
//        event("vtk_file");
//    }
//
///**
//We set a constant viscosity corresponding to a Reynolds number of 40, 100,
//based on the cylinder diameter (1) and the inflow velocity (1). */
//
//    event properties (i++){
//    frame_BP (fs, 0);
//    foreach_face() muv.x[] = mu1;
//    boundary((scalar *){muv});
//    update_targetU(target_Uv, t+dt);
//    }
//
//    event set_dtmax (i++){
//    NITERMIN = 5;
//    NITERMAX = 100;
//    DT *= 1.05;
//    DT = min(DT, 0.01);
//    fprintf(ferr, "set_dtmax: tnext= %g Dt= %g", tnext, DT);
//    }
//
//    event advection_term(i = 0){
////    u.n[right] = neumann(0);
////    u.t[right] = neumann(0);
////    uold.n[right] = neumann(0);
////    uold.t[right] = neumann(0);
//        u.n[right] = outxBC(x,y,t,dt);
//        u.t[right] = outyBC(x,y,t,dt);
//        uold.n[right] = outxBC(x,y,t,dt);
//        uold.t[right] = outyBC(x,y,t,dt);
//    }
//
//    event advection_term(i++){
////    uf.n[right]  = outxBC(x,y,t,0.5*dt);
////    uf.t[right]  = outyBC(x,y,t,0.5*dt);
////    uf.n[right]  = neumann(0);
////    uf.t[right]  = neumann(0);
//    //deep copy of uold with ghost cells
//    foreach() foreach_dimension() uold.x[] = u.x[];
//    for (int b = 0; b < nboundary; b++) foreach_boundary(b) foreach_dimension() uold.x[ghost] = u.x[ghost];
//    boundary({uold});
//    }
//
//    event projection(i++){
////    uf.n[right]  = outxBC(x,y,t,dt);
////    uf.t[right]  = outyBC(x,y,t,dt);
////    uf.n[right]  = neumann(0);
////    uf.t[right]  = neumann(0);
//    }
//
////event end_timestep(i++){
////    double xref = -1.5, yref = -1.5, zref = 0; //reference location
////    double pref = 0;                                            //reference pressure
////    double pcor = interpolate (p, xref, yref, zref) - pref;
////    fprintf(ferr, "pressure correction=%g at x=%g y=%g z=%g\n", pcor, xref, yref, zref);
////#if _MPI
////    MPI_Allreduce (MPI_IN_PLACE, &pcor, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
////#endif
////    foreach() p[] -= pcor;
////    boundary ({p});
////}
//
//    event end_timestep(i++){
//    double Cx=0, Cy=0;
//    foreach_boundary(left, reduction(+:Cx) reduction(+:Cy)){
//    Cx += uf.x[]*Delta;
//    Cy += uf.y[]*Delta;
//    }
//    C.x = Cx/Ldomain;
//    C.y = Cy/Ldomain;
//    fprintf(ferr, "mean mass velocity at left: C=%g %g\n", C.x, C.y);
//
//    Cx=0, Cy=0;
//    foreach_boundary(right, reduction(+:Cx) reduction(+:Cy)){
//    Cx += uf.x[ghost]*Delta;
//    Cy += uf.y[ghost]*Delta;
//    }
//    fprintf(ferr, "mean mass velocity at right: C=%g %g\n", Cx/Ldomain, Cy/Ldomain);
//    }
//    event end_timestep (i += 1){
////event end_timestep (i += 100){
//    vector uexact[];
//    scalar pexact[];
//    vorticity (u, omega);
//    double Luinf = 0, Lpinf = 0, du, dp, maxp = 0, maxu = 0;
//    int curlevel=1;
//    theory(uexact, pexact, t+dt, RE);
//    scalar l[]; foreach(reduction(max:curlevel)) {l[] = level; if (curlevel<level) curlevel=level;}
//#ifdef DEBUG_MODE
//    vector gradu[], gradv[];
//    gradients({u.x, u.y}, {gradu, gradv});
//    foreach(){
//        mod_du_dx[] = sqrt(sq(gradu.x[]) + sq(gradu.y[]) + sq(gradv.x[]) + sq(gradv.y[]));
//    }
//#endif
//    foreach(reduction(max:Luinf) reduction(max:Lpinf) reduction(max:maxu) reduction(max:maxp)){
//    if (fs[] == 0){ //in inner frame
//    du = sqrt(sq(u.x[] - uexact.x[]) + sq(u.y[] - uexact.y[]));
//    dp = fabs(p[] - pexact[]);
//    if (du > Luinf) Luinf = du;
//    if (dp > Lpinf) Lpinf = dp;
//    if (norm(uexact) > maxu) maxu = norm(uexact);
//    if (fabs(pexact[]) > maxp) maxp = fabs(pexact[]);
//    }
//    }
//
//    int tnc = count_cells();
//    int nmax = (int)pow (2, maxlevel*dimension);
//    fprintf (ferr, "i= %d t+dt= %g dt= %g Luinf= %g Lpinf= %g Luinf_rel= %g Lpinf_rel= %g eta_s= %g count_cells= %d nmax= %d compress_ratio= %g curlevel= %d\n", i, t+dt, dt, Luinf, Lpinf, Luinf/u0, Lpinf/p0, eta_s, tnc, nmax, (double)tnc/nmax, curlevel);
//    double eps_arr[] = {1, 1, 1, 1, 1, 1, 1};
//    MinMaxValues((scalar *){p, pexact, u.x, u.y, uexact.x, uexact.y, omega}, eps_arr);
//    }
//
////Output
//    event vtk_file (i++){
////event vtk_file (i+=10){
////event vtk_file (t += 1){
//    char subname[80]; sprintf(subname, "couette_not_periodic_adapt_jac2");
//    scalar l[], pid_num[]; foreach() {l[] = level; pid_num[] = pid();}
//    vector uexact[];
//    scalar pexact[];
//    theory(uexact, pexact, t+mydt, RE);
//    vector uf_low[], uf_up[], u_low[], u_up[], p_low[], p_up[];
//    face vector uf_exact[];
//
//    foreach(){
//        foreach_dimension()
//        {
//            uf_low.x[] = uf.x[];
//            uf_up.x[]  = uf.x[1];
//        }
//        u_low.x[] = u.x[-1];
//        u_low.y[] = u.y[-1];
//        u_up.x[] = u.x[1];
//        u_up.y[] = u.y[1];
//        p_low.x[] = p[-1,0];
//        p_low.y[] = p[0,-1];
//        p_up.x[] = p[1,0];
//        p_up.y[] = p[0,1];
//    }
//    boundary((scalar *){uf_low, uf_up, u_low, u_up, p_low, p_up});
//#ifdef DEBUG_MODE
//    double d = 0;
//    foreach() {
//        d = 0.;
//        foreach_dimension(){
//            d += uf.x[1] - uf.x[];
//        }
//        divutmpAfter[] = d/(dt*Delta);
//    }
//    boundary((scalar*){divutmpAfter});
//#endif
//    output_vtu_MPI( (scalar *) {phi, fs, omega, p, pf, pexact, l, divutmpAfter, pid_num, deltap, divutmp},
//    (vector *) {u, uexact, uf_low, uf_up, deltag, g, target_Uv, u_low, u_up, p_low, p_up, dbp, total_rhs, tmp_err_u, divtauu}, subname, t + mydt );
////    output_vtu_MPI( (scalar *) {fs, omega, p, pexact, l, pid_num}, (vector *) {u, uexact}, subname, periodic_BC);
//    fprintf(ferr, "dt=%g mydt=%g\n", dt, mydt);
//    }
///**
//We adapt according to the error on the embedded geometry, velocity*/
////#define ADAPT_SCALARS {omega, mod_du_dx}
////#define ADAPT_EPS_SCALARS (double[]){eps*2*pi, eps*sqrt(2.0)*pi}
////event adapt (i++){
//////    double eps_arr[] = ADAPT_EPS_SCALARS;//??? sometimes doesn't work
//////    MinMaxValues(ADAPT_SCALARS, eps_arr);
//////    adapt_wavelet ((scalar *) {fs, omega}, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
////    adapt_wavelet ((scalar *) {u.x, u.y}, (double[]){eps, eps}, maxlevel = maxlevel, minlevel = minlevel);
////}
//
////event stop (i = 10){
//    event stop (t = 50){
//        event("end_timestep");
//    };
//
///**
//## See also
//
//* [Same example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/cylinder.html)
//*/
