#define BRINKMAN_PENALIZATION 1
#define DEBUG_BRINKMAN_PENALIZATION 1
#define DEBUG_OUTPUT_VTU_MPI
//#define FILTERED
//#define JACOBI 1
#define EPS_MAXA 2
scalar omega[];
face vector fs_face[];
#include "../src_local/centered-weugene.h"
#include "view.h"
#include "../src_local/output_vtu_foreach.h"
#include "../src_local/three-phase-weugene.h"
#include "tension.h"
#include "fracface.h"

int maxlevel = 10;
int minlevel = 4;
//double U0=0.01, rhol=1e+3, sig=73e-3, Lchar=5e-3, mul=1e-3, Lb=0.3, Rb=0.125;
//double Rrho=1000, Rmu=53.73;
double U0=1, rhol=1, sig=0.0005, Lchar=1, mul=1, Lb=0.3, Rb=0.0625;
double Rrho=1, Rmu=1;
double RE, CA, Rrad;
double xs0 = -0.4, rad = 0.0625;
coord vc = {1.0, 0.0, 0.0};
double signvc = 1;


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
    eta_s = 1e-5;
    TOLERANCE = 1e-4;
    NITERMAX = 30;
    N = 1 << minlevel;
    periodic(top);
    periodic(right);

    Rb = Rrad*rad;
    rho1 = 1.; rho2 = rho1/Rrho; rho3 = max(rho1, rho2);
    mu1 = 1./RE; mu2 = mu1/Rmu; mu3 = max(mu1, mu2);
    f.sigma = 1./RE/CA;
    signvc = (vc.x > 0) ? 1 : (vc.x < 0)? -1 : 0;
    fprintf(ferr, "maxlevel=%d tol=%g NITERMAX=%d\n"
                  "RE=%g CA=%g Rb/rad=%g rho1/rho2=%g mu1/mu2=%g\n"
                  "mu1=%g mu2=%g rho1=%g rho2=%g sigma=%g\n",
            maxlevel, TOLERANCE, NITERMAX,
            RE, CA, Rrad, Rrho, Rmu,
            mu1, mu2, rho1, rho2, f.sigma);

    const vector U_sol[] = {vc.x, vc.y, vc.z};
    target_U = U_sol;
    run();
}

scalar divu[];
#define tmpposx (xs0 + vc.x*t)
#define outOfBox ((positionX < X0) || (positionX > X0 +L0))
double posx (double t){
    double positionX = tmpposx;
    while(outOfBox){
        positionX -= signvc*L0;
    }
    return positionX;
}
void soild_fs(scalar fs, double t){
    double deltaT=0.0;
    double tt = max(t-deltaT,0);
//    vc.x = sin(t);
//    const vector U_sol[] = {vc.x, vc.y, vc.z};
//    target_U = U_sol;
//    foreach() {fprintf(ferr, "targetU=%g %g vc.x=%g\n", target_U.x[], target_U.y[], vc.x);break;}
    vertex scalar phi[];
    face vector ff[];
    double x00 = posx(tt);
    foreach_vertex() {
        phi[] = HUGE;
        for (int xi=-L0; xi <=L0; xi +=L0) {
            double x1 = x00 + xi;
            phi[] = intersection(phi[], (sq(x - x1) + sq(y - vc.y * tt) - sq(rad)));
        }
        phi[] = -phi[];
    }
    boundary ({phi});
    fractions (phi, fs, ff);
    fs.refine = fs.prolongation = fraction_refine;
    boundary({fs});
}
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
            soild_fs (fs, 0);
            bubbles(f);
            foreach() {
                u.x[] = U0*fs[];
                u.y[] = 0;
            }
            boundary({f,fs,u});
            foreach_face() fs_face.x[] = 0.5*(fs[-1] + fs[]);
            //call viscosity
        }while (adapt_wavelet({f, fs, u}, (double []){1e-5, 1e-5, 1e-2, 1e-2}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
        refine((sq(x-xs0) + sq(y) <sq(1.2*rad)) && (sq(x-xs0) + sq(y) >sq(0.8*rad)) && level <10);
        DT=1e-9;
    }
}

event set_dtmax (i++) if (i<500) DT *= 1.05;

event vof(i++){
//    foreach_face() uf.x[] *= 1 - fs_face.x[];
}


//event moving_cylinder (i++) {
//    soild_fs(fs, t);
//}
event properties (i++) {
    soild_fs(fs, t + dt);
//    face_fraction (fs, fs_face);
    foreach_face() fs_face.x[] = 0.5*(fs[-1] + fs[]);
}

//event viscous_term (i++){
//    soild_fs(fs, t + dt);
//}
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
    char subname[80]; sprintf(subname, "rk");
    scalar l[];
    foreach() {l[] = level;}

    vector mapped_data_lower[], mapped_data_upper[], mu_low[];
    vector fs_lower[], fs_upper[];
    foreach() {
        foreach_dimension()
        {
            mapped_data_lower.x[] = uf.x[];
            mapped_data_upper.x[] = uf.x[1];
            fs_lower.x[] = fs_face.x[];
            fs_upper.x[] = fs_face.x[1];
            mu_low.x[]=mu.x[];
        }
    }
    output_vtu_MPI( (scalar *) {fs, f, omega, p, l, divu, rho}, (vector *) {u, g, dbp, total_rhs, a, mu_low, mapped_data_lower, mapped_data_upper, fs_lower, fs_upper}, subname, 1 );
}

#define ADAPT_SCALARS {f, fs, omega}
#define ADAPT_EPS_SCALARS {1e-5, 1e-4, 1e-2}
event adapt (i++){
    double eps_arr[] = ADAPT_EPS_SCALARS;
    MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
    fs.refine = fs.prolongation = fraction_refine;
    boundary({fs});
}

event snapshot (t += 0.1) {
      char name[80];
      sprintf (name, "snapshot-%g", t);
      dump (name);
}

event stop(t = 10);










//#define BRINKMAN_PENALIZATION 1
//#define DEBUG_BRINKMAN_PENALIZATION 1
//#define DEBUG_OUTPUT_VTU_MPI
////#define FILTERED
////#define JACOBI
//#include "../src_local/centered-weugene.h"
//#include "view.h"
//#include "../src_local/output_vtu_foreach.h"
//#include "two-phase.h"
//#include "tension.h"
//
//int maxlevel = 11;
//int minlevel = 4;
//double U0=1, rhol=1, sig=0.01, Lchar=1, mul=0.00078125, Lb=0.5, Rb=0.1;
//double RE, CA;
//double Rrho=1, Rmu=1;
//double xs0 = 0, rad = 0.0625;
//coord vc = {1, 0.0, 0.0};
//scalar fs[], omega[];
//
///**
//The domain is the periodic unit square centered on the origin. */
//
///**
//The fluid is injected on the left boundary with a unit velocity. The
//tracer is injected in the lower-half of the left boundary. An outflow
//condition is used on the right boundary. */
//
////u.n[left]  = dirichlet(0);
////u.t[left]  = dirichlet(0);
////p[left]    = neumann(0);
////pf[left]   = neumann(0);
////
////u.n[right] = neumann(0.);
////u.t[right] = neumann(0.);
////p[right]   = neumann(0.);
////pf[right]  = neumann(0.);
//
//int main(int argc, char * argv[])
//{
//    if (argc > 1) {
//        maxlevel = atoi(argv[2]); //convert from string to float
//    }
//    size (8.0);
//    origin (-0.5, -L0/2.);
//    eta_s = 1e-10;
//    DT = 1e-2;
//    TOLERANCE = 1e-8;
//    NITERMAX = 150;
//    N = 512;
//    periodic(top);
//    periodic(right);
//    RE=U0*(2*rad)*rhol/mul; CA=U0*mul/sig;
//    rho1 = 1.; rho2 = rho1/Rrho;
//    mu1 = 1./RE; mu2 = mu1/Rmu;
//    f.sigma = 1./RE/CA;
//    fprintf(ferr, "RE=%g CA=%g \n"
//                  "mu1=%g mu2=%g rho1=%g rho2=%g sigma=%g\n",
//            RE, CA, mu1, mu2, rho1, rho2, f.sigma);
//
//    const vector U_wall[] = {vc.x, vc.y, vc.z};
//    target_U = U_wall;
//    run();
//}
//
//scalar divu[];
//void soild_fs(scalar fs, double t){
//    fraction (fs, sq(rad) - sq(x - xs0 - vc.x*t) - sq(y - vc.y*t));
//    fs.refine = fs.prolongation = fraction_refine;
//    boundary({fs});
//}
//event init (t = 0) {
//    if (!restore (file = "restart")) {
//        int it = 0;
//        do {
//            fraction(f, (sq(x - Lb - xs0) + sq(y) > sq(Rb))? 1 : -1);
//            soild_fs (fs, 0);
//            foreach() u.x[] = U0*fs[];
//            boundary(all);
//        }while (adapt_wavelet({f, fs}, (double []){1e-3, 1e-3}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
//    }
//}
//
//event moving_cylinder (i++) {
//soild_fs(fs, t);
//}
//
//void correct_press(scalar p, int i){
//    double press = 0;
//    int ip = 0;
//#if 1 // Left bottom Corner
//    foreach(){
//        if (ip == 0){
//            press = p[];
//            ip++;
//            @if _MPI
//                MPI_Bcast(&press, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//            @endif
//        }
//    }
//#else //average value
//    press = normf(p).avg;
//#endif
//    foreach(){
//        p[] -= press;
//    }
//}
//
//event end_timestep(i++){
//correct_press(p, i);
//}
//
//event logfile (i++) {
//foreach() {
//    divu[] = 0;
//    foreach_dimension() divu[] += (uf.x[1]-uf.x[])/Delta;
//}
//double Linfu = -10;
//foreach( reduction(max:Linfu) ){
//if (fabs(divu[]) > Linfu) Linfu = fabs(divu[]);
//}
//fprintf (ferr, "i=%d t=%g dt=%g iter_p=%d iter_u=%d div u=%g \n", i, t, dt, mgp.i, mgu.i, Linfu);
//}
//
///**
//We produce animations of the vorticity and tracer fields... */
//
//event images (t += 0.1) {
//static FILE * fp = popen ("ppm2gif > vort.gif", "w");
//vorticity (u, omega);
///**
//Cells for which *m* is negative will be black in the movie. */
//scalar m[];
//foreach()
//        m[] = 0.5 - fs[];
//                boundary ({m});
//output_ppm (omega, fp, box = {{-0.5,-0.5},{7.5,0.5}}, mask = m,
//        min=-10, max=10, linear=true);
//}
//
////Output
////event vtk_file (i += 1){
//event vtk_file (t += 0.01){
//char subname[80]; sprintf(subname, "rk");
//scalar l[];
//vorticity (u, omega);
//foreach() {l[] = level; omega[] *= 1 - fs[];}
//output_vtu_MPI( (scalar *) {fs, f, omega, p, l, divu}, (vector *) {u, dbp, total_rhs}, subname, 1 );
//}
//
//#define ADAPT_SCALARS {f, fs, u}
//#define ADAPT_EPS_SCALARS {1e-4, 1e-4, 3e-2, 3e-2}
//event adapt (i++){
//double eps_arr[] = ADAPT_EPS_SCALARS;
////	MinMaxValues(ADAPT_SCALARS, eps_arr);
//adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
//fs.refine = fs.prolongation = fraction_refine;
//boundary({fs});
//}
//
//event stop(t = 7);
