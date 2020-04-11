#define BRINKMAN_PENALIZATION 1
#define DEBUG_BRINKMAN_PENALIZATION 1
#define DEBUG_OUTPUT_VTU_MPI
//#define FILTERED
#include "../src_local/centered-weugene.h"
#include "view.h"
#include "../src_local/output_vtu_foreach.h"
#include "two-phase.h"
#include "tension.h"

int maxlevel = 8;
int minlevel = 4;
double U0 = -1, rhol=1, sig=0.01, mul=0.00078125, Lb=0.3, Rb=0.125;
double RE, CA;
double Rrho=1, Rmu=1;
double xx0 = -0.1, rad = 0.0625;
scalar fs[], omega[], divu[];

/**
The domain is the periodic unit square centered on the origin. */

/**
The fluid is injected on the left boundary with a unit velocity. The
tracer is injected in the lower-half of the left boundary. An outflow
condition is used on the right boundary. */

//u.n[left]  = neumann(0.);
//p[left]    = dirichlet(0.);
//pf[left]   = dirichlet(0.);
//
//
//u.n[right] = dirichlet(U0);
//p[right]   = neumann(0.);
//pf[right]  = neumann(0.);
void soild_fs(scalar fs, double t){
    fraction (fs, sq(rad) - sq(x - xx0) - sq(y));
    fs.refine = fs.prolongation = fraction_refine;
    boundary({fs});
}

int main(int argc, char * argv[])
{
    if (argc > 1) {
        maxlevel = atoi(argv[2]); //convert from string to float
    }
    size (1.0);
    origin (-0.5*L0, -0.5*L0);
    eta_s = 1e-10;
    DT = 1e-2;
    TOLERANCE = 1e-8;
    NITERMAX = 150;
    N = 512;
    periodic(top);
    periodic(right);
    RE=fabs(U0)*(2*rad)*rhol/mul; CA=fabs(U0)*mul/sig;
    rho1 = 1.; rho2 = rho1/Rrho;
    mu1 = 1./RE; mu2 = mu1/Rmu;
    f.sigma = 1./RE/CA;
    fprintf(ferr, "RE=%g CA=%g \n"
                  "mu1=%g mu2=%g rho1=%g rho2=%g sigma=%g\n",
            RE, CA, mu1, mu2, rho1, rho2, f.sigma);
    run();
}

scalar un[];

event init (t = 0) {
    if (!restore (file = "restart")) {
        int it = 0;
        do {
            fraction(f, (sq(x - Lb - xx0) + sq(y) > sq(Rb))? 1 : -1);
            soild_fs (fs, 0);
            foreach() u.x[] = U0*(1 - fs[]);
            boundary(all);
        }while (adapt_wavelet({fs, f}, (double []){1e-4, 1e-4}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
        boundary(all);
    }
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
}

event logfile (i++) {
foreach() {
    divu[] = 0;
    foreach_dimension() divu[] += (uf.x[1]-uf.x[])/Delta;
}
double Linfu = -10;
foreach( reduction(max:Linfu) ){
if (fabs(divu[]) > Linfu) Linfu = fabs(divu[]);
}
fprintf (ferr, "i=%d t=%g dt=%g iter_p=%d iter_u=%d div u=%g \n", i, t, dt, mgp.i, mgu.i, Linfu);
}

/**
We produce animations of the vorticity and tracer fields... */

event images (t += 0.1) {
    static FILE * fp = popen ("ppm2gif > vort.gif", "w");
    vorticity (u, omega);
    output_ppm (omega, fp, min=-10, max=10, linear=true);
}

//Output
//event vtk_file (i += 1){
event vtk_file (t += 0.01){
    char subname[80]; sprintf(subname, "rk");
    scalar l[];
    vorticity (u, omega);
    foreach() {l[] = level;}

    #if DEBUG_BRINKMAN_PENALIZATION!=1
        output_vtu_MPI( (scalar *) {fs, f, omega, p, l, divu}, (vector *) {u}, subname, 1 );
    #else
        output_vtu_MPI( (scalar *) {fs, f, omega, p, l, divu}, (vector *) {u, dbp, total_rhs}, subname, 1 );
    #endif
}

#define ADAPT_SCALARS {fs, f, u}
#define ADAPT_EPS_SCALARS {1e-4, 1e-4, 3e-2, 3e-2}
event adapt (i++){
    double eps_arr[] = ADAPT_EPS_SCALARS;
    //	MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
    fs.refine = fs.prolongation = fraction_refine;
    boundary({fs});
}
event stop(t = 3);







//#define BRINKMAN_PENALIZATION 1
//#define DEBUG_BRINKMAN_PENALIZATION 1
//#define DEBUG_OUTPUT_VTU_MPI
////#define FILTERED
//#include "../src_local/centered-weugene.h"
//#include "view.h"
//#include "../src_local/output_vtu_foreach.h"
//#include "two-phase.h"
//#include "tension.h"
//
//int maxlevel = 11;
//int minlevel = 4;
//double U0 = -1, rhol=1, sig=0.01, Lchar=1, mul=0.00078125, Lb=0.5, Rb=0.1;
//double RE, CA;
//double Rrho=1, Rmu=1;
//double xx0 = 6, rad = 0.0625;
//scalar fs[], omega[], divu[];
//
///**
//The domain is the periodic unit square centered on the origin. */
//
///**
//The fluid is injected on the left boundary with a unit velocity. The
//tracer is injected in the lower-half of the left boundary. An outflow
//condition is used on the right boundary. */
//
////u.n[left]  = neumann(0.);
////p[left]    = dirichlet(0.);
////pf[left]   = dirichlet(0.);
////
////
////u.n[right] = dirichlet(U0);
////p[right]   = neumann(0.);
////pf[right]  = neumann(0.);
//void soild_fs(scalar fs, double t){
//    fraction (fs, sq(rad) - sq(x - xx0) - sq(y));
//    fs.refine = fs.prolongation = fraction_refine;
//    boundary({fs});
//}
//
//int main(int argc, char * argv[])
//{
//	if (argc > 1) {
//		maxlevel = atoi(argv[2]); //convert from string to float
//	}
//	size (8.0);
//	origin (-0.5, -L0/2.);
//	eta_s = 1e-10;
//    DT = 1e-2;
//	TOLERANCE = 1e-8;
//    NITERMAX = 150;
//	N = 512;
//	periodic(top);
//	periodic(right);
//    RE=fabs(U0)*(2*rad)*rhol/mul; CA=fabs(U0)*mul/sig;
//    rho1 = 1.; rho2 = rho1/Rrho;
//    mu1 = 1./RE; mu2 = mu1/Rmu;
//    f.sigma = 1./RE/CA;
//    fprintf(ferr, "RE=%g CA=%g \n"
//                  "mu1=%g mu2=%g rho1=%g rho2=%g sigma=%g\n",
//            RE, CA, mu1, mu2, rho1, rho2, f.sigma);
//	run();
//}
//
//scalar un[];
//
//event init (t = 0) {
//	if (!restore (file = "restart")) {
//		int it = 0;
//		do {
//            fraction(f, (sq(x - Lb - xx0) + sq(y) > sq(Rb))? 1 : -1);
//            soild_fs (fs, 0);
//            foreach() u.x[] = U0*(1 - fs[]);
//            boundary(all);
//		}while (adapt_wavelet({fs, f}, (double []){1e-3, 1e-3}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
//	    boundary(all);
//	}
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
//    correct_press(p, i);
//}
//
//event logfile (i++) {
//    foreach() {
//        divu[] = 0;
//        foreach_dimension() divu[] += (uf.x[1]-uf.x[])/Delta;
//    }
//    double Linfu = -10;
//    foreach( reduction(max:Linfu) ){
//        if (fabs(divu[]) > Linfu) Linfu = fabs(divu[]);
//    }
//    fprintf (ferr, "i=%d t=%g dt=%g iter_p=%d iter_u=%d div u=%g \n", i, t, dt, mgp.i, mgu.i, Linfu);
//}
//
///**
//We produce animations of the vorticity and tracer fields... */
//
//event images (t += 0.1) {
//    static FILE * fp = popen ("ppm2gif > vort.gif", "w");
//    vorticity (u, omega);
//    /**
//    Cells for which *m* is negative will be black in the movie. */
//    scalar m[];
//    foreach()
//            m[] = 0.5 - fs[];
//                    boundary ({m});
//    output_ppm (omega, fp, box = {{-0.5,-0.5},{7.5,0.5}}, mask = m,
//            min=-10, max=10, linear=true);
//}
//
////Output
////event vtk_file (i += 1){
//event vtk_file (t += 0.01){
//    char subname[80]; sprintf(subname, "rk");
//    scalar l[];
//    vorticity (u, omega);
//    foreach() {l[] = level; omega[] *= 1 - fs[]; }
//
//    #if DEBUG_BRINKMAN_PENALIZATION!=1
//        output_vtu_MPI( (scalar *) {fs, f, omega, p, l, divu}, (vector *) {u}, subname, 1 );
//    #else
//        output_vtu_MPI( (scalar *) {fs, f, omega, p, l, divu}, (vector *) {u, dbp, total_rhs}, subname, 1 );
//    #endif
//}
//
//#define ADAPT_SCALARS {fs, f, u}
//#define ADAPT_EPS_SCALARS {1e-4, 1e-4, 3e-2, 3e-2}
//event adapt (i++){
//    double eps_arr[] = ADAPT_EPS_SCALARS;
//    //	MinMaxValues(ADAPT_SCALARS, eps_arr);
//    adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
//    fs.refine = fs.prolongation = fraction_refine;
//    boundary({fs});
//}
//event stop(t = 7);