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
double mydt = 0;
int maxlevel;
int minlevel;
double Ldomain = 4.0, RE = 30., rho1 = 1, G = 1, h = 1, U = 1;
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
        phi[] = ( y <= 1.001 && y>=-0.001) ? 0 : 1;
    }
    boundary ({phi});
    fractions (phi, fs, fs_face);
}

#define uxe(x,y,t) (((G/(2*mu.x[]))*y*(h-y) + U*y/h)*(1.0 - fs[]) + fs[]*target_Uv.x[]);
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
u.t[left]  = neumann(0); // available x,y,z at the computational domain
u.n[left]  = neumann(0);
p[left]    = dirichlet(pe(x,y,t));
//pf[left]   = dirichlet(pe(x,y,t));
//p[left]    = dirichlet(pe(x,y,t+dt));
//pf[left]   = dirichlet(pe(x,y,t+0.5*dt));

//u.t[right] = dirichlet(uye(x,y,t+dt));
//u.n[right] = dirichlet(uxe(x,y,t+dt));
u.t[right] = neumann(0);
u.n[right] = neumann(0);
p[right]   = dirichlet(pe(x,y,t));
//pf[right]  = dirichlet(pe(x,y,t));
//p[right]   = dirichlet(pe(x,y,t+dt));
//pf[right]  = dirichlet(pe(x,y,t+0.5*dt));

u.t[bottom] = dirichlet(0);
u.n[bottom] = dirichlet(0);
p[bottom] = neumann(0);
//pf[bottom]= neumann(0);
//p[bottom] = dirichlet(pe(x,y,t+dt));
//pf[bottom]= dirichlet(pe(x,y,t+0.5*dt));

u.t[top] = dirichlet(1);
u.n[top] = dirichlet(0);
p[top]   = neumann(0);
//pf[top]  = neumann(0);
//p[top]   = dirichlet(pe(x,y,t+dt));
//pf[top]  = dirichlet(pe(x,y,t+0.5*dt));

//uf.n[left]   = neumann(0);
//uf.t[left]   = neumann(0);
//uf.n[right]  = neumann(0);
//uf.t[right]  = neumann(0);
uf.n[bottom] = dirichlet(0);
uf.t[bottom] = dirichlet(0);
uf.n[top]    = dirichlet(0);
uf.t[top]    = dirichlet(1);

deltap[left]   = dirichlet(0);
deltap[right]  = dirichlet(0);
deltap[bottom] = neumann(0);
deltap[top]    = neumann(0);

int main(int argc, char * argv[]) {
    minlevel = 4;
    maxlevel = 10;
    eta_s=1e-4;
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
    DT = 1e-2;
    CFL = 0.1;
    TOLERANCE = 1e-8;
    RELATIVE_RES_TOLERANCE = 0.1;
    NITERMAX = 30;
    mu = muv;
    stokes = true;

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
            foreach() u.x[] = fs[]*target_Uv.x[];
            boundary({u});
        } while ( ++it <= 10 && adapt_wavelet((scalar *){fs, u}, (double []){eps, eps, eps}, maxlevel=maxlevel, minlevel=minlevel).nf != 0);
        foreach() {
            p[] = 0;//pexact[];// my initial guess
            pf[] = 0;//pexact[];// my initial guess
        }
        boundary({p, pf});
    }
    event("vtk_file");

}

/**
We set a constant viscosity corresponding to a Reynolds number of 40, 100,
based on the cylinder diameter (1) and the inflow velocity (1). */

event properties (i++)
{
    frame_BP (fs, 0);
    foreach_face() muv.x[] = 1.0/RE;
    boundary((scalar *){muv});
    update_targetU(target_Uv, t+dt);
}

event set_dtmax (i++) {
    NITERMIN=1;
    NITERMAX=100;

    DT *= 1.05;
//    DT = min(DT, CFL*Ldomain/pow(2, maxlevel+3));

    fprintf(ferr, "set_dtmax: tnext= %g Dt= %g", tnext, DT);
}


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
/**
We produce animations of the vorticity and tracer fields... */

//event movies (i += 4; t <= 15.)
//{
//    scalar omega[], m[];
//    vorticity (u, omega);
//    foreach()
//            m[] = fs[] - 0.5;
//    boundary ({m});
//    output_ppm (omega, file = "vort.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
//        min = -10, max = 10, linear = true, mask = m);
//#if TURN_ON_TRACER == 1
//    output_ppm (f, file = "f.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
//        linear = false, min = 0, max = 1, mask = m);
//#endif
//}

//Output
//event vtk_file (i++){
event vtk_file (i+=5){
//event vtk_file (t += 0.01){
    char subname[80]; sprintf(subname, "vortex_BP_parallel_not_periodic_cc");
    scalar l[], pid_num[]; foreach() {l[] = level; pid_num[] = pid();}
    vector uexact[];
    scalar pexact[];
    theory(uexact, pexact, t+mydt, RE);
    vector uf_low[], uf_up[], uf_exact_low[], uf_exact_up[];
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

//    foreach_boundary(left){
//        fprintf(ferr, "++left:%g %g %g %g du=%g x=%g %g\n ", uf_exact.x[-1], uf_exact.x[], uf_exact.x[1], uxe(x,y,t+mydt), uf_exact.x[]-uxe(x,y,t+mydt), x,y);
//        assert (uf_exact.x[] == uxe(x,y,t+mydt));
//    }
//
//    foreach_boundary(right){
//        fprintf(ferr, "++right:%g %g %g %g x=%g %g\n ", uf_exact.x[-1], uf_exact.x[], uf_exact.x[1], uxe(x,y,t+mydt), x,y);
//        assert (uf_exact.x[ghost] == uxe(x,y,t+mydt));
//    }
//    foreach_boundary(bottom){
//        fprintf(ferr, "++bottom:%g %g %g %g x=%g %g\n ", uf_exact.y[0,-1], uf_exact.y[], uf_exact.y[0,1], uye(x,y,t+mydt), x,y);
//        assert (uf_exact.y[] == uye(x,y,t+mydt));
//    }
//    foreach_boundary(top){
//        fprintf(ferr, "++top:%g %g %g %g x=%g %g\n ", uf_exact.y[0,-1], uf_exact.y[], uf_exact.y[0,1], uye(x,y,t+mydt), x,y);
//        assert (uf_exact.y[ghost] == uye(x,y,t+mydt));
//    }
    foreach(){
        foreach_dimension()
        {
            uf_low.x[] = uf.x[];
            uf_up.x[]  = uf.x[1];
            uf_exact_low.x[] = uf_exact.x[];
            uf_exact_up.x[]  = uf_exact.x[1];
        }
    }
    boundary((scalar *){uf_low, uf_up, uf_exact_low, uf_exact_up});
    output_vtu_MPI( (scalar *) {fs, omega, p, pf, pexact, l, divutmpAfter, pid_num, deltap, divutmp},
            (vector *) {u, uexact, uf_low, uf_up, uf_exact_low, uf_exact_up, deltag, g}, subname, periodic_BC );
//    output_vtu_MPI( (scalar *) {fs, omega, p, pexact, l, pid_num}, (vector *) {u, uexact}, subname, periodic_BC);
    fprintf(ferr, "dt=%g mydt=%g\n", dt, mydt);

//    foreach_boundary(left){
////            fprintf(ferr, "left:%g %g %g x=%g %g\n ", uf.x[], uf_exact.x[], uxe(x,y,t+mydt), x,y);
//            assert (uf.x[] == uf_exact.x[]);
//    }
//    foreach_boundary(right){
////            fprintf(ferr, "right:%g %g %g x=%g %g\n ", uf.x[1], uf_exact.x[1], uxe(x,y,t+mydt), x,y);
//            assert (uf.x[ghost] == uf_exact.x[ghost]);
//    }
//    foreach_boundary(bottom){
////        fprintf(ferr, "bottom:%g %g %g x=%g %g\n ", uf.y[], uf_exact.y[], uye(x,y,t+mydt), x,y);
//        assert (uf.y[] == uf_exact.y[]);
//    }
//    foreach_boundary(top){
////        fprintf(ferr, "top:%12.8g %12.8g %12.8g %12.8g %12.8g %12.8g %12.8g x=%g %g\n ",  uf.y[-1], uf.y[], uf.y[1], uf_exact.y[-1],  uf_exact.y[], uf_exact.y[1], uye(x,y,t+mydt), x,y);
//        assert (uf.y[ghost] == uf_exact.y[ghost]);
//    }
//    return 1;
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
event stop (t = 300){
    event("end_timestep");
};
//event stop (i = 10);
/**
## See also

* [Same example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/cylinder.html)
*/
