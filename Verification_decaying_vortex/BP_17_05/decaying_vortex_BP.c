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
#define DEBUG_OUTPUT_VTU_MPI
#define DEBUG_MINMAXVALUES
#define FILTERED
#define JACOBI 1
#define EPS_MAXA 2
#define RELATIVE_RESIDUAL
#define MODIFIED_CHORIN 1
scalar omega[], fs[];
face vector fs_face[];
vector target_Uv[];
face vector target_Ufv[];
scalar divutmp[];
face vector my_u_rhs[];
face vector my_alpha[];
double eta_chorin = 1e-6;
#include "../src_local/centered-weugene.h"
//#include "view.h"
#include "../src_local/output_vtu_foreach.h"
#include "fractions.h"
face vector muv[];
vector uexact[];
scalar pexact[];
int maxlevel = 8;
int minlevel = 4;
int dtlimiter = 0;
double Ldomain = 3.0, RE = 30.;
double ueps = 1e-3, fseps = 1e-5, peps=1e-3;

void frame_BP(scalar fs, face vector face_vector, double t){
    vertex scalar phi[];
    foreach_vertex() {
        phi[] = ( fabs(x) <= 1 && fabs(y) <= 1 ) ? -1 : 1;
    }
    boundary ({phi});
    fractions (phi, fs, face_vector);
    foreach_face(){
        fs_face.x[] = face_value(fs, 0);
    }
    boundary((scalar *){fs, fs_face});
}

#define uxe (-cos(pi*x) * sin(pi*y) * exp(-2.0 * sq(pi) * t / RE))
#define uye ( sin(pi*x) * cos(pi*y) * exp(-2.0 * sq(pi) * t / RE))
#define pxe (-0.25 *  (cos(2*pi*x) + cos(2.0*pi*y)) * exp(-4.0 * sq(pi) * t / RE))
void theory(vector u, scalar p, double t, double RE){
    foreach(){
        u.x[] = uxe;
        u.y[] = uye;
        p[]   = pxe;
    }
    boundary({u, p});
    fprintf(ferr, "theory  t= %g\n", t);
}
//$$u(x, y, t) = − \cos 􏱲x \sin 􏱲y \exp{−2􏱲\pi^2 t/Re}$$
//$$v(x, y, t) =   \sin 􏱲x \cos 􏱲y \exp{−2􏱲\pi^2 t/Re}$$
//$$p(x, y, t) =  −\frac14(\cos 2􏱲x + \cos 2􏱲y)\exp{−4\pi^2 t/Re}$$
//in 2D each cells must have 2 BC, in 3D - 3BC
u.t[left]  = neumann(0);
p[left]    = dirichlet(pxe);
pf[left]   = dirichlet(pxe);

u.t[right] = neumann(0);
p[right]   = dirichlet(pxe);
pf[right]  = dirichlet(pxe);

u.t[top] = neumann(0);
p[top]   = dirichlet(pxe);
pf[top]  = dirichlet(pxe);

u.t[bottom] = neumann(0);
p[bottom]   = dirichlet(pxe);
pf[bottom]  = dirichlet(pxe);


int main(int argc, char * argv[]) {
    eta_s = 1e-6;
    if (argc > 1) {
        maxlevel = atoi(argv[1]); //convert from string to float
    }
    if (argc > 2) {
        eta_s = atof(argv[2]); //convert from string to float
    }
    if (argc > 3) {
        dtlimiter = atoi(argv[3]); //convert from string to float
    }

    size (Ldomain);
    origin (-0.5*Ldomain, -0.5*Ldomain);
    DT = 1e-8;
    CFL = 0.4;
    TOLERANCE = 1e-8;
    RELATIVE_RES_TOLERANCE = 0.1;
    NITERMAX = 30;

    mu = muv;

    target_U = target_Uv;
    target_Uf = target_Ufv;
//    for(maxlevel=7; maxlevel<=11; maxlevel++) {
        N = 1<<maxlevel;
        fprintf(ferr, "maxlevel=%d DTlimiter=%d eta=%g\n", maxlevel, dtlimiter, eta_s);
        fprintf(ferr, "TOL=%g NITERMAX=%d Re=%g Rel_res_tol=%g\n", TOLERANCE, NITERMAX, RE, RELATIVE_RES_TOLERANCE);
        fprintf(ferr, "CFL=%g feps=%g peps=%g ueps=%g\n", CFL, fseps, peps, ueps);
        run();
//    }
}

void update_targetU(vector target_Uv, face vector target_Uf, double t){
    foreach() {
        target_Uv.x[] = uxe;
        target_Uv.y[] = uye;
    }

    foreach_face(x) target_Uf.x[] = uxe;
    foreach_face(y) target_Uf.y[] = uye;
    boundary((scalar *){target_Uv, target_Uf});
}

event init (t = 0)
{
    if (!restore (file = "restart")) {
        int it = 0;
        do {
            frame_BP (fs, fs_face, 0);
            theory(uexact, pexact, 0, RE);
        } while ( ++it <= 10 && adapt_wavelet((scalar *){fs, pexact, uexact}, (double []){fseps, peps, ueps, ueps}, maxlevel=maxlevel, minlevel=minlevel).nf != 0);

        foreach() {
            p[] = pexact[];// useless
            foreach_dimension() {u.x[] = uexact.x[];}
        }
        boundary({p, u});
        update_targetU(target_Uv, target_Ufv, 0);
    }

}

/**
We set a constant viscosity corresponding to a Reynolds number of 40, 100,
based on the cylinder diameter (1) and the inflow velocity (1). */

event properties (i++)
{
    frame_BP (fs, fs_face, 0);
    foreach_face() muv.x[] = 1.0/RE;
    boundary((scalar *){muv});
    update_targetU(target_Uv, target_Ufv, t+dt);
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
    if(dtlimiter) DT = min(DT, eta_s);
    fprintf(ferr, "set_dtmax: tnext= %g Dt= %g", tnext, DT);
}

/**
We check the number of iterations of the Poisson and viscous
problems. */
scalar ptmp[];
void correct_press(scalar p, scalar fs, int i){
    double press = 0;
    int ip = 0;
#if 0 // Left bottom Corner
    foreach(){
        if (ip == 0){
            press = p[];
            ip++;
            break;
        }
    }
#else //average value
    foreach(reduction(+:press)){
        if (fabs(fs[]) < SEPS) {
            press += p[]*dv()*(1. - fs[]);
        }
    }
    @if _MPI
        MPI_Bcast(&press, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    @endif
    press /=4; // inner frame are is 4
#endif

    foreach(){
        p[] -= press;
    }
//    boundary((scalar *){p});
    fprintf(ferr, "correct_press= %g \n", press);
}


event end_timestep (i++){
//    correct_press(p, fs, i);
    vorticity (u, omega);
    double Luinf = 0, Lpinf = 0, du, dp, maxp = 0, maxu = 0;
    theory(uexact, pexact, t+dt, RE);
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
    fprintf (ferr, "i= %d t+dt= %g dt= %g Luinf= %g Lpinf= %g Luinf_rel= %g Lpinf_rel= %g iter_p= %d iter_u= %d \n", i, t+dt, dt, Luinf, Lpinf, Luinf/maxu, Lpinf/maxp, mgp.i, mgu.i);
    double eps_arr[] = {1, 1, 1, 1, 1, 1, 1, 1};
    MinMaxValues((scalar *){p, pexact, u.x, u.y, uexact.x, uexact.y, omega, divutmp}, eps_arr);
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
event vtk_file (t += 0.01){
    char subname[80]; sprintf(subname, "vortex_BP");
    scalar l[], omega[];
    foreach() {l[] = level;}
    vector my_u_rhs_low[], my_u_rhs_up[], my_alpha_low[], my_alpha_up[], fs_face_low[], fs_face_up[];
    foreach(){
        foreach_dimension()
        {
            my_u_rhs_low.x[] = my_u_rhs.x[];
            my_u_rhs_up.x[] = my_u_rhs.x[1];
            my_alpha_low.x[] = my_alpha.x[];
            my_alpha_up.x[] = my_alpha.x[1];
            fs_face_low.x[] = fs_face.x[];
            fs_face_up.x[] = fs_face.x[1];
        }
    }
    boundary((scalar*){my_u_rhs_low, my_u_rhs_up, my_alpha_low, my_alpha_up});
    output_vtu_MPI( (scalar *) {fs, omega, p, pexact, l, divutmp},
            (vector *) {u, uexact, my_u_rhs_low, my_u_rhs_up, my_alpha_low, my_alpha_up, fs_face_low, fs_face_up}, subname, 0 );
}
/**
We adapt according to the error on the embedded geometry, velocity*/
#define ADAPT_SCALARS {fs, p, u}
#define ADAPT_EPS_SCALARS {fseps, peps, ueps, ueps}
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
