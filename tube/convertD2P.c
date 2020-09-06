//#include "run.h"
//#include "timestep.h"
scalar fs[], omega[], l2[];
#include "grid/octree.h"

#include "run.h"
#include "timestep.h"
#include "bcg.h"
#if EMBED
#include "viscosity-embed.h"
#else
#ifndef BRINKMAN_PENALIZATION
#include "viscosity.h"
#else
#include "../src_local/viscosity-weugene.h"
#include "../src_local/utils-weugene.h"
#endif
#endif

//#include "../src_local/centered-weugene.h"
//#include "two-phase.h"
//#include "tension.h"
#include "utils.h"
#include "lambda2.h"
#include "../src_local/output_vtu_foreach.h"
#include "../src_local/adapt_wavelet_limited.h"
#include "../src_local/adapt2.h"
#include "maxruntime.h"
#include <stdio.h>
#include <stdlib.h>
#include <wordexp.h>
#include <ctype.h>
double Unean = 1;

double fseps = 1e-3, ueps = 1e-2;
double length_min = 1e+30, length_max = -1e+30, length = 1;

//scalar p[], pf[];
vector u[];
//vector g[];
//face vector uf[];

(const) face vector mu = zerof, a = zerof, alpha = unityf, kappa = zerof;
(const) scalar rho = unity;
mgstats mgp, mgpf, mgu;
double dtmax;

scalar f[], * interfaces = {f};
double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;
//face vector alphav[];
//scalar rhov[];


int maxlevel = 10;
int minlevel = 5;
int LEVEL = 6;
int i_take = 1;
int adapt_method = 0;
double lDomain;
double myt=0;
wordexp_t fp;
char **w, dump_name[30];

int maXlevel(double x,double y, double z){
    double x0 = fabs(x - 0.5*(length_min + length_max));
    int n = ceil(max(0, 0.5*(x0/length - 3)));
    return max(maxlevel-n, 10);
}
//get time from dump name (dump-1.1234 => 1.1234)
double get_double(const char *str)
{
    /* First skip non-digit characters */
    /* Special case to handle negative numbers and the `+` sign */
    while (*str && !(isdigit(*str) || ((*str == '-' || *str == '+') && isdigit(*(str + 1)))))
        str++;
    /* The parse to a double */
    return strtod(str, NULL);
}

int main (int argc, char * argv[]) {
    maxruntime (&argc, argv);
    if (argc > 1)
        maxlevel = atoi (argv[1]);
    if (argc > 2)
        lDomain = atof (argv[2]); // set L0
    if (argc > 3)
        i_take = atoi (argv[3]); // set which dump files will be converted: each $(i_take)th

    size(lDomain);
    origin(0., -L0/2., -L0/2.);
    init_grid(1 << LEVEL);

    wordexp("dump-*", &fp, 0);
    w = fp.we_wordv;
    for (int i = 0; i < fp.we_wordc; i++) fprintf(ferr, "All dump files: %s\n", w[i]);

    for (int i = 0; i < fp.we_wordc; i += i_take) {
//        dt = 0;
        myt =  fabs(get_double(w[i]));
        strcpy(dump_name, w[i]);
        fprintf(ferr, "reading dump file: %s at time t= %g\n", dump_name, myt);
        run();
    }
    wordfree(&fp);
}

event init (t = 0) {
    bool success = restore (file = dump_name);
    fprintf(ferr, "file has been read: L0=%g, maxlevel=%d\n", L0, maxlevel);
    if (L0 != lDomain) {
        fprintf(ferr, "L0=%g doesn't coinside with input lDomain=%g\n", L0, lDomain);
        return 0;
    }
    if (!success) {
        fprintf(ferr, "can't open the file %s. Missing this file, go to the next file\n", dump_name);
        return 0;
    }
}

event coarsen_grid(i++){
    double xcg = 0, dvtmp, volume = 0, volumeg = 0 ;
    length_min = 1e+30, length_max = -1e+30, length = 0;
    foreach( reduction(+:xcg) reduction(+:volume) reduction(+:volumeg)) {
        if (fs[]<1){
            dvtmp = (1.0 - f[])*(1.0 - fs[])*dv(); // gas volume
            volumeg += dvtmp;//gas liquid
            volume += (1.0 - fs[])*dv();//channel volume
            xcg   += x*dvtmp;// Along x
        }
    }
    xcg /= volumeg;
    length_min = xcg - 5;
    length_max = xcg + 4;
    length = length_max - length_min;


    fprintf (ferr, "x= %g length_min= %g length_max= %g length= %g it_fp= %d\n"
                   "volume= %g volumeg= %g\n",
             xcg, length_min, length_max, length, iter_fp,
             volume, volumeg);
}



event vtk_file (i++)
{
    char subname[80]; sprintf(subname, "dump2pvd_compressed");
    scalar l[]; foreach() l[] = level;
    vorticity (u, omega);
    lambda2 (u, l2);
//    unrefine ( x < 15 && level >= 1);
    unrefine ( (x < length_min || x > length_max) && level >= 1);
    unrefine ( (sq(y) + sq(z) > sq(0.55)) && level >= 1);
//    unrefine ( (sq(y) + sq(z) > sq(0.51)) && level >= 1);
//    unrefine ( (x < length_min || x > length_max) && (sq(y) + sq(z) > sq(0.51)) && level >= 1);
    output_vtu_MPI( subname, myt, (scalar *) {fs, f, l, l2, omega}, (vector *) {u});
}

//#define ADAPT_SCALARS {f, u.x, u.y, u.z}
//#define ADAPT_SCALARS_SMOOTHED {f, u.x, u.y, u.z}
//#define ADAPT_EPS_SCALARS {fseps, ueps, ueps, ueps}
//#define ADAPT_MAXLEVEL {maxlevel, max(maxlevel-2,10), max(maxlevel-2,10), max(maxlevel-2,10)}
//event adapt(i++){
//    if (adapt_method == 0)
//        adapt_wavelet ((scalar *) ADAPT_SCALARS, (double []) ADAPT_EPS_SCALARS, maxlevel = maxlevel, minlevel = minlevel);
//    else if (adapt_method == 1)
//    //        adapt_wavelet_limited  ((scalar *) {f, u_mag}, (double []) {fseps, ueps}, maXlevel, minlevel);
//        adapt_wavelet_limited  ((scalar *) ADAPT_SCALARS, (double []) ADAPT_EPS_SCALARS, maXlevel, minlevel);
//    else if (adapt_method == 2)
//        adapt_wavelet2((scalar *)ADAPT_SCALARS, (double []) ADAPT_EPS_SCALARS,(int []){maxlevel, maxlevel-1, maxlevel-2, maxlevel-2, maxlevel-2}, minlevel);
//
//    return 0;
//}
event stop(t = 100);






//#define BRINKMAN_PENALIZATION 1
//#define DEBUG_MINMAXVALUES
////#define DEBUG_BRINKMAN_PENALIZATION
//#define DEBUG_MODE_POISSON
////#define DEBUG_OUTPUT_VTU_MPI
//#define FILTERED
//#define JACOBI 1
////#define PRINT_ALL_VALUES
//#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
////#define STOKES
//double myt=0;
//scalar fs[];
//scalar omega[];
//scalar l2[];
//#include "grid/octree.h"
//#include "../src_local/centered-weugene.h"
//#include "two-phase.h"
//#ifdef STOKES
//#include "navier-stokes/conserving.h"
//#endif
//#include "tension.h"
//#include "../src_local/adapt_wavelet_limited.h"
//#include "../src_local/adapt2.h"
//#include "../src_local/utils-weugene.h"
//#include "utils.h"
//#include "lambda2.h"
//#include "../src_local/output_vtu_foreach.h"
//#include "maxruntime.h"
//#include <stdio.h>
//#include <stdlib.h>
//#include <wordexp.h>
//#include <ctype.h>
//
//#define AIR_WATER
////#define AIR_GLYCEROL
//#define uexact(x,y,z) 2.*(1. - 4*sq(y) - 4*sq(z))
////#define uexact(x,y,z) 0.25*(G/mu1)*(sq(0.5) - sq(y) - sq(z))
////Channel cross section Lyy*Lzz
//double Vd, Vdst, deq, dst = 0.2, rst = 0.1, r_bub, l_bub;
//double RhoR, MuR;
//#if defined(AIR_WATER)
//double Rho1 = 997, Rho2 = 1.204;
//double Mu1 = 0.88e-3, Mu2 = 0.019e-3;
//double Sigma = 72.8e-3;
//double diam_tube = 514e-6;
//double dt_vtk = 1e-2;
//double lDomain = 20;
//#elif defined(AIR_GLYCEROL)
//double Rho1 = 1250, Rho2 = 1.204;
//    double Mu1 = 550e-3, Mu2 = 0.019e-3;
//    double Sigma = 63.4e-3;
//    double diam_tube = 494e-6;
//    double dt_vtk = 1e-3;
//    double lDomain = 10;
//#endif
//double Ca; // Ca = Mu*Ud/sigma
//double Ca_mod; // Ca_mod = Mu*Umean/sigma
//double Re; //Reynolds
//double G;
//double Umean;
//double x_init = 2;
//int maxlevel = 10;
//int minlevel = 5;
//int LEVEL = 6;
//int adapt_method = 1; // 0 - traditional, 1 - using limitation, 2 - using array for maxlevel
//int i_take = 1;
//double fseps = 1e-3, ueps = 5e-2;
//bool ellipse_shape = false, cylinder_shape = false;
//scalar un[];
////char * dump_name = "dump-0000";
//int main (int argc, char * argv[]) {
//    maxruntime (&argc, argv);
//
//    eta_s = 1e-5;
//    TOLERANCE = 1e-6;
//    NITERMIN = 1;
//    NITERMAX = 100;
//    DT = 1e-3;
//    cylinder_shape = true;
////    resolve_capillary_effects = true;
////    relative_residual_poisson = true;
////    relative_residual_viscous = true;
//    fs.refine = fs.prolongation = fraction_refine;
//// Case 9e
//    Ca = 0.163; //Ca = Ud*Mu1/sigma
//    Umean = 0.01145; // m/s
//    Vd = 0.0780e-9; // m^3
//
//// Case 10e
////    Ca = 0.023; //Ca = Ud*Mu1/sigma
////    Umean = 1.580; // m/s
////    Vd = 0.2179e-9; // m^3
//
//    if (argc > 1)
//        maxlevel = atoi (argv[1]);
//    if (argc > 2)
//        lDomain = atof (argv[2]);
//    if (argc > 3)
//        i_take = atoi (argv[3]);
//    if (argc > 4)
//        iter_fp = atoi (argv[4]);
//
//    size(lDomain);
//    origin(0., -L0/2., -L0/2.);
//    init_grid(1 << LEVEL);
//
//    deq = pow(6*Vd/pi, 1./3.);// 0.0005301091821 m
//    dst = deq/diam_tube;// 1.0730955104
//    rst = 0.5*dst;
//    Vdst = (4./3.)*pi*cube(rst);
////    Umean = G*sq(0.5)/(8*Mu1);
//    Ca_mod = Mu1*Umean/Sigma;
//    Re = Umean*diam_tube*Rho1/Mu1;
//    G = 32.0*Mu1*Umean/sq(diam_tube);
//
//    if (ellipse_shape || dst < 0.9) {
//        r_bub = min(rst, 0.4);
//        l_bub = cube(rst) / sq(r_bub);
//        ellipse_shape = true;
//        cylinder_shape = false;
//    }else if (cylinder_shape){
//        r_bub = 0.45;
//        l_bub = (Vdst - (4./3.)*pi*cube(r_bub))/(pi*sq(r_bub));
//    } else{
//        assert(false && "set shape");
//    }
//    x_init = 1.7*l_bub;
//
//    fprintf(ferr,"BP:             eta_s=%g,     DT=%g\n"
//                 "Solver:         NITERMIN=%d   NITERMAX=%d      TOLERANCE=%g  relative_residual_poisson=%d relative_residual_viscous=%d\n"
//                 "OUTPUT:         dt_vtk=%g i_take=%d iter_fp=%d\n"
//                 "ADAPT:          minlevel=%d,  maxlevel=%d      adapt_meth=%d fseps=%g ueps=%g\n"
//                 "Properties(SI): Mu1=%g Mu2=%g Rho1=%g Rho2=%g  Sigma=%g G=%g Umean=%g\n"
//                 "Apparatus:      diam_tube=%g  tube_length=%g\n"
//                 "Bubble:         Vd=%g deq=%g  ellipse_shape=%d cylinder_shape=%d\n",
//            eta_s, DT,
//            NITERMIN, NITERMAX, TOLERANCE, relative_residual_poisson, relative_residual_viscous,
//            dt_vtk, i_take, iter_fp,
//            minlevel, maxlevel, adapt_method, fseps, ueps,
//            Mu1, Mu2, Rho1, Rho2, Sigma, G, Umean,
//            diam_tube, L0,
//            Vd, deq, ellipse_shape, cylinder_shape );
//    // Dimensionless parameters:
//    // Averaging on diam_tube=1 and Umean=1, Mu1=1 and Rho1=1 p' = p/(Rho1*Umean^2)
////    G /= Mu1*Umean/sq(diam_tube);
//    G /= Rho1*sq(Umean)/diam_tube;
//    Umean /= Umean;
//    diam_tube /= diam_tube;
//    RhoR = Rho1/Rho2;
//    MuR = Mu1/Mu2;
//
//    rho1 = 1.;// water
//    rho2 = 1./RhoR; // air
//    mu1 = 1./Re;
//    mu2 = 1./(MuR*Re);
//    f.sigma = 1./(Re*Ca_mod);
//    fprintf(ferr,"Dimensionless Parameters: mu1=%g mu2=%g rho1=%g rho2=%g sigma=%g G=%g  Umean=%g\n"
//                 "Dimensionless nums:       Re=%g  Ca=%g  Ca_mod=%g\n"
//                 "Bubble:                   Vdst=%g dst=%g  rst=%g  r_bub=%g l_bub=%g x_init=%g\n",
//            mu1, mu2, rho1, rho2, f.sigma, G, Umean,
//            Re, Ca, Ca_mod,
//            Vdst, dst, rst, r_bub, l_bub, x_init);
//    run();
//}
//
//
//double get_double(const char *str)
//{
//    /* First skip non-digit characters */
//    /* Special case to handle negative numbers and the `+` sign */
//    while (*str && !(isdigit(*str) || ((*str == '-' || *str == '+') && isdigit(*(str + 1)))))
//        str++;
//
//    /* The parse to a double */
//    return strtod(str, NULL);
//}
//
//event init (t = 0) {
//    wordexp_t fp;
//    char **w;
//
//    wordexp("dump-*", &fp, 0);
//    w = fp.we_wordv;
//    for (int i = 0; i < fp.we_wordc; i++) fprintf(ferr, "dump file: %s\n", w[i]);
//    for (int i = 0; i < fp.we_wordc; i += i_take) {
//        dt = 0;
//        myt =  fabs(get_double(w[i]));
//        fprintf(ferr, "reading dump file: %s at time t= %g\n", w[i], myt);
//        bool success = restore (file = w[i]);
//        fprintf(ferr, "file is read: L0=%g, maxlevel=%d\n", L0, maxlevel);
//        if (L0 != lDomain) {
//            fprintf(ferr, "L0=%g doesn't coinside with input lDomain=%g\n", L0, lDomain);
//            continue;
//        }
//        if (!success) {
//            fprintf(ferr, "can't open the file %s. Missing this file, go to the next file\n", w[i]);
//            continue;
//        }else{
//            fprintf(ferr, "dump2pvtu ... %s\n", w[i]);
//            event("vtk_file");
//        }
//    }
//    wordfree(&fp);
//    exit(EXIT_SUCCESS);
//    exit(1);
//}
//
//event coarsen_grid(i+=10000){
//double avggas = L0*1 - normf(f).avg;
//scalar umag[];
//double xcg = 0, volume = 0, volumeg = 0, velamean = 0, velgx = 0, velgy = 0, dvtmp, min_r=0, mean_r=0, max_r=0, r=0;
//double length_min = 0, length_max = 0, length = 0;
//foreach(reduction(+:xcg) reduction(+:volume) reduction(+:volumeg)
//reduction(+:velgx) reduction(+:velgy) reduction(+:velamean)
//reduction(min:min_r) reduction(max:max_r) ) {
//if (fs[]<1){
//dvtmp = (1.0 - f[])*(1.0 - fs[])*dv(); // gas volume
//volumeg += dvtmp;//gas liquid
//volume += (1.0 - fs[])*dv();//channel volume
//umag[] = norm(u); // the length of u
//xcg   += x*dvtmp;// Along x
//velamean += umag[]*dvtmp;//mean velocity of gas
//velgx += u.x[]*dvtmp;//mean velocity of gas Ox
//velgy += u.y[]*dvtmp;//mean velocity of gas Oy
//if (f[] < 1 ) {
//r = sqrt(sq(y) + sq(z)); // accurate order of Delta/2
//if (r < min_r) min_r = r;
//if (r > max_r) max_r = r;
//length = x; // accurate order of Delta/2
//if (length < length_min) length_min = length;
//if (length > length_max) length_max = length;
//}
//}
//}
//mean_r = 0.5*(min_r + max_r);
//xcg /= volumeg; velgx /= volumeg; velgy /= volumeg; velamean /= volumeg;
////norm statu = normf_weugene(umag, fs); // outputs avg, rms, max, volume
//fprintf (ferr, "maxlevel= %d i= %d t= %g dt= %g avggas= %g velgx= %g valgy= %g velamean= %g velgx/U0-1= %g \n"
//"length_min= %g length_max= %g xcg= %g thickness_min= %g thickness_mean= %g thickness_max= %g length= %g it_fp= %d\n",
//maxlevel, i, t, dt, avggas, velgx, velgy, velamean, (velgx/Umean - 1),
//length_min, length_max, xcg, 0.5 - max_r, 0.5 - mean_r, 0.5 - min_r, length_max - length_min, iter_fp);
////    refine ( ((x > length_min - 4) && (x < length_max + 4) && (sq(y) + sq(z) < sq(0.51))) && level <= 5);
////    unrefine ( sq(y) + sq(z) > sq(0.51)  && level >= 1);
////    unrefine (level >= 5);
//}
////event vtk_file (i += 1)
//event vtk_file (t += dt_vtk)
//{
//char subname[80]; sprintf(subname, "dump2pvd_compressed");
//scalar l[]; foreach() l[] = level;
//        vorticity (u, omega);
//lambda2 (u, l2);
//event("coarsen_grid");
//fprintf(ferr, "haha\n");
//
//output_vtu_MPI( subname, myt, (scalar *) {p, fs, f, l, residual_of_p, l2, omega}, (vector *) {u});
//}
//
//event adapt(i++){
//unrefine(x > 2 || x < -2);
//return 0;
//}
//event stop(t=L0/Umean);