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
//#define DEBUG_MINMAXVALUES
#define DEBUG_MODE
#define FILTERED
//#define JACOBI 1
#define EPS_MAXA 2
#define RELATIVE_RESIDUAL
#define MODIFIED_CHORIN 0
//#define PERIODIC_BC
scalar omega[], fs[];
vector target_Uv[];
vector deltag[];
scalar deltap[];
scalar phi[];
vector uold[];
vector tmp_err_u[];
tensor myt = new tensor;
double mydt = 0;
int maxlevel;
int minlevel;
double Ldomain = 4.0, RE = 30., rho1 = 1, G = 1, h = 1, U = 1, C, mu1;
double eps = 1e-2;
#ifdef DEBUG_MODE
    scalar divutmp[], divutmpAfter[];
    scalar mod_du_dx[];
    vector conv_term[];
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
//        phi[] = 0;
//        phi[] = -( pow(x,10) + pow(y,10) - 0.117 );
//        phi[] = ( fabs(x) <= 1 && fabs(y) <= 1 ) ? 0 : 1;
        phi[] = ( y*(y-1) + Ldomain/pow(2, maxlevel+1));
    }
    boundary ({phi});
    fractions (phi, fs, fs_face);
}

//max U = 4.25
#define uxe(x,y,t) (((G/(2.0*mu1))*y*(h-y) + U*y/h)*(1.0 - fs[]) + fs[]*target_Uv.x[])
#define uye(x,y,t) 0
#define pe(x,y,t)  G*(0.5*Ldomain - x)
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
u.n[left]   = uxe(x,y,t); // available x,y,z at the computational domain
u.n[right]  = neumann(0);
u.n[bottom] = dirichlet(0);
u.n[top]    = dirichlet(0);

//u.t[left]   = dirichlet(0);
u.t[left]   = 0;
u.t[right]  = neumann(0);
u.t[bottom] = dirichlet(0);
u.t[top]    = dirichlet(1);

//uold.n[left]   = dirichlet(uxe(x,y,t)); // available x,y,z at the computational domain
uold.n[left]   = uxe(x,y,t); // available x,y,z at the computational domain
uold.n[right]  = neumann(0);
uold.n[bottom] = dirichlet(0);
uold.n[top]    = dirichlet(0);

//uold.t[left]   = dirichlet(0);
uold.t[left]   = 0;
uold.t[right]  = neumann(0);
uold.t[bottom] = dirichlet(0);
uold.t[top]    = dirichlet(1);

p[left]    = neumann(G);
p[right]   = dirichlet(pe(x,y,t));
p[bottom]  = dirichlet(pe(x,y,t));
p[top]     = dirichlet(pe(x,y,t));

//p[left]    = neumann(-rho1*(uxe(x,y,t+dt) - uf.x[])/dt);
//p[right]   = neumann(-rho1*(uxe(x,y,t+dt) - uf.x[ghost])/dt);
//p[bottom] = neumann(-rho1*(uye(x,y,t+dt) - uf.y[])/dt);
//p[top] = neumann(-rho1*(uye(x,y,t+dt) - uf.y[])/dt);

uf.n[left]   = uxe(x,y,t);
uf.n[right]  = neumann(0);
uf.n[bottom] = 0;
uf.n[top]    = 0;

// du/dt = F(x,y,t,C=1,u^n)
uf.t[left]   = 0;
uf.t[right]  = neumann(0);
uf.t[bottom] = 0;
uf.t[top]    = 1;
//u.t[right] = neumann(0);
//u.n[right] = neumann(0);
//u.n[right] = neumann(-u.x[ghost]/(C*dt))/(1-Delta/(C*dt));
//p[right]   = neumann(-G);

//uf.n[right]  = neumann(0);
//uf.t[right]  = neumann(0);
//uf.n[right]  = dirichlet(uxe(x,y,t+dt));
//uf.t[right]  = dirichlet(uxe(x,y,t+dt));

deltap[left]   = neumann(0);
deltap[right]  = dirichlet(0);
deltap[bottom] = neumann(0);
deltap[top]    = neumann(0);

phi[left]   = neumann(uf.x[] - uxe(x,y,t));
phi[right]  = neumann(uf.x[] - uxe(x,y,t));
phi[bottom] = neumann(uf.y[] - uye(x,y,t));//?
phi[top]    = neumann(uf.y[] - uye(x,y,t));

int main(int argc, char * argv[]) {
    minlevel = 6;
    maxlevel = 6;
    eta_s=1e-3;
//    cells_per_zone = 2;
    eps = 1e-3;
    if (argc > 1) {
        maxlevel = atoi(argv[1]);
    }
    if (argc > 2) {
        eta_s = atof(argv[2]);
    }
    if (argc > 3) {
        eps = atof(argv[3]);
    }
#ifdef PERIODIC_BC
    periodic(top);
#endif
    size (Ldomain);
    origin (-0.5*Ldomain, -0.5*Ldomain);
    DT = 1e-3;
    CFL = 0.1;
    TOLERANCE = 1e-10;
    RELATIVE_RES_TOLERANCE = 0.1;
    NITERMAX = 30;
    mu1 = 1./RE;
    mu = muv;
    stokes = true;
    C = (G*sq(h) + 6*U*mu1)/(12.*mu1);
    target_U = target_Uv;
    N = 1<<minlevel;
    fprintf(ferr, "maxlevel= %d eta= %g eps= %g\n", maxlevel, eta_s, eps);
    fprintf(ferr, "TOL=%g NITERMAX=%d Re=%g Rel_res_tol=%g\n", TOLERANCE, NITERMAX, RE, RELATIVE_RES_TOLERANCE);
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
        do {
            event ("properties");
            theory(uexact, pexact, 0, RE);
            foreach() foreach_dimension() u.x[] = uexact.x[];
//            foreach() u.x[] = fs[]*target_Uv.x[];
            boundary({u});
        } while ( ++it <= 10 && adapt_wavelet((scalar *){fs, u}, (double []){eps, eps, eps}, maxlevel=maxlevel, minlevel=minlevel).nf != 0);
        foreach() {
            p[] = pexact[];// my initial guess
            pf[] = 0;//pexact[];// my initial guess
            foreach_dimension() uold.x[] = u.x[];
        }
        boundary({p, pf, uold});
    }
    event("vtk_file");

}

/**
We set a constant viscosity corresponding to a Reynolds number of 40, 100,
based on the cylinder diameter (1) and the inflow velocity (1). */

event properties (i++)
{
    frame_BP (fs, 0);
    foreach_face() muv.x[] = mu1;
    boundary((scalar *){muv});
    update_targetU(target_Uv, t+dt);
}

event set_dtmax (i++) {
    NITERMIN=1;
    NITERMAX=100;
    DT *= 1.05;
    DT = min(DT, 0.001);
    fprintf(ferr, "set_dtmax: tnext= %g Dt= %g", tnext, DT);
}

event advection_term(i++){
    u.n[left]    = uxe(x,y,t+dt);
    uold.n[left] = uxe(x,y,t+dt);
    uf.n[left]   = uxe(x,y,t+dt);
//    u.n[right]  = (uold.x[] - u.x[] + u.x[ghost] + 2.0*(dt*C/Delta))/(1.0 + 2.0*dt*C/Delta);
//    uold.n[right]  = (uold.x[] - u.x[] + u.x[ghost] + 2.0*(dt*C/Delta))/(1.0 + 2.0*dt*C/Delta);
    foreach() foreach_dimension() uold.x[] = u.x[];
    boundary({uold});
}
event end_timestep (i += 100){
    vector uexact[];
    scalar pexact[];
    vorticity (u, omega);
    double Luinf = 0, Lpinf = 0, du, dp, maxp = 0, maxu = 0, d;
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
#ifdef DEBUG_MODE
    foreach() {
        d = 0.;
        foreach_dimension(){
            d += uf.x[1] - uf.x[];
        }
        divutmpAfter[] = d/(dt*Delta);
    }
    boundary((scalar*){divutmpAfter});
#endif
    int tnc = count_cells();
    int nmax = (int)pow (2, maxlevel*dimension);
    fprintf (ferr, "i= %d t+dt= %g dt= %g Luinf= %g Lpinf= %g Luinf_rel= %g Lpinf_rel= %g eta_s= %g count_cells= %d nmax= %d compress_ratio= %g curlevel= %d\n", i, t+dt, dt, Luinf, Lpinf, Luinf/u0, Lpinf/p0, eta_s, tnc, nmax, (double)tnc/nmax, curlevel);
    double eps_arr[] = {1, 1, 1, 1, 1, 1, 1};
    MinMaxValues((scalar *){p, pexact, u.x, u.y, uexact.x, uexact.y, omega}, eps_arr);
}

//Output
event vtk_file (i++){
//event vtk_file (i+=5){
//event vtk_file (t += 1){
    char subname[80]; sprintf(subname, "coette");
    scalar l[], pid_num[]; foreach() {l[] = level; pid_num[] = pid();}
    vector uexact[];
    scalar pexact[];
    theory(uexact, pexact, t+mydt, RE);
    vector uf_low[], uf_up[], uf_exact_low[], uf_exact_up[], u_low[], u_up[], p_low[], p_up[];
    face vector uf_exact[];

    foreach_face(x){
        uf_exact.x[] = uxe(x,y,t+mydt);
    }
    foreach_face(y){
        uf_exact.y[] = uye(x,y,t+mydt);
    }
    uf_exact.t[left]   = uye(x,y,t+mydt); // face location at x
    uf_exact.n[left]   = uxe(x,y,t+mydt);
    uf_exact.t[right]  = uye(x,y,t+mydt);
    uf_exact.n[right]  = uxe(x,y,t+mydt);
    uf_exact.t[bottom] = uxe(x,y,t+mydt);
    uf_exact.n[bottom] = uye(x,y,t+mydt);
    uf_exact.t[top]    = uxe(x,y,t+mydt);
    uf_exact.n[top]    = uye(x,y,t+mydt);
    boundary((scalar *){uf_exact});

    foreach(){
        foreach_dimension()
        {
            uf_low.x[] = uf.x[];
            uf_up.x[]  = uf.x[1];
            uf_exact_low.x[] = uf_exact.x[];
            uf_exact_up.x[]  = uf_exact.x[1];
        }
        u_low.x[] = u.x[-1];
        u_low.y[] = u.y[-1];
        u_up.x[] = u.x[1];
        u_up.y[] = u.y[1];
        p_low.x[] = p[-1,0];
        p_low.y[] = p[0,-1];
        p_up.x[] = p[1,0];
        p_up.y[] = p[0,1];
    }
    boundary((scalar *){uf_low, uf_up, uf_exact_low, uf_exact_up, u_low, u_up, p_low, p_up});
    output_vtu_MPI( (scalar *) {phi, fs, omega, p, pf, pexact, l, divutmpAfter, pid_num, deltap, divutmp},
            (vector *) {u, uexact, uf_low, uf_up, uf_exact_low, uf_exact_up, deltag, g, target_Uv, u_low, u_up, p_low, p_up, dbp, total_rhs, tmp_err_u, muv}, subname, periodic_BC );
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

event stop (i = 5){
//event stop (t = 100){
    event("end_timestep");
};
//event stop (i = 10);
/**
## See also

* [Same example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/cylinder.html)
*/
