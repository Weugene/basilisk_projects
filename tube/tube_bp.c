#define BRINKMAN_PENALIZATION 1
#define DEBUG_MINMAXVALUES
//#define DEBUG_BRINKMAN_PENALIZATION
//#define DEBUG_MODE_POISSON
//#define DEBUG_OUTPUT_VTU_MPI
#define FILTERED
#define JACOBI 1

//#define PRINT_ALL_VALUES
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
//#define STOKES

scalar fs[];
face vector fs_face=zerof;
#include "grid/octree.h"
#include "centered-weugene.h"
#include "two-phase.h"
#ifdef STOKES
    #include "navier-stokes/conserving.h"
#endif
#include "tension.h"
#include "adapt_wavelet_limited.h"
#include "adapt2.h"
#include "utils-weugene.h"
#include "utils.h"
#include "lambda2.h"
#if 0
    #include "output_vtu_foreach.h"
#else
    int iter_fp = 0;
#endif
#include "maxruntime.h"
#include "tag.h"

#define uexact(x,y,z) 2.*(1. - 4*sq(y) - 4*sq(z))
//#define uexact(x,y,z) 0.25*(G/mu1)*(sq(0.5) - sq(y) - sq(z))
//Channel cross section Lyy*Lzz
double Vd, Vdst, deq, dst = 0.2, rst = 0.1, r_bub, l_bub;
double RhoR, MuR;
double Rho1, Rho2;
double Mu1, Mu2;
double Sigma;
double diam_tube;
double dt_vtk;
double lDomain;
int bubcase;

typedef struct {
    double Ca;    // Ca=mu*Ub/sigma
    double Vd;    // volume
    double Uc;    // Umean
    double Ud;    // Ububble
    double delta; // thickness
} Cases;
#define zv {0,0,0,0,0}
//         CA       V_d[l]    Uc[m/s]  Ud[m/s] delta*[-]
Cases cases[28]={
        zv, // 0 is empty
        zv,zv,zv,zv,zv,zv,zv, // 1-7 Air-glycerol
        {6.44e-4,  0.0349e-9, 0.0454, 0.0533, 0.1055},// 8  Air-Water delta?
        {8.42e-4,  0.0439e-9, 0.0593, 0.0697, 0.074}, // 9  Air-Water delta?
        {7.34e-4,  0.0602e-9, 0.0557, 0.0607, 0.027}, // 10 Air-Water delta?
        {6.64e-4,  0.0745e-9, 0.0534, 0.0550, 0.001}, // 11 Air-Water delta?
        {6.697e-4, 0.1893e-9, 0.0543, 0.0554, 0.006}, // 12 Air-Water
        zv,zv,zv,zv,zv, // 13-17 Air-glycerol
        {0.003, 0.1751e-9, 0.242,  0.261, 0.013},  // 18 Air-Water
        {0.008, 0.1715e-9, 0.666,  0.704, 0.023},  // 19 Air-Water
        {0.0098, 0.2208e-9, 0.757, 0.815, 0.025}, // 20 Air-Water
        {0.015,  0.1882e-9, 1.118 , 1.293, 0.039},  // 21 Air-Water
        {0.024,  0.2179e-9, 1.580, 1.944, 0.054},  // 22 Air-Water
        {0.034,  0.2179e-9, 2.060, 2.511, 1e-9},  // 23 Air-Water Ud, delta - garbage
        {0.0455, 0.2179e-9, 2.575, 3.165, 1e-9},  // 24 Air-Water Ud, delta - garbage
        {0.056,  0.2179e-9, 3.09,  4.65,   1e-9},  // 25 Air-Water Ud, delta - garbage
        {0.065,  0.2179e-9, 3.602, 5.40,   1e-9},  // 26 Air-Water Ud, delta - garbage
        {0.074,  0.2179e-9, 4.117, 6.17,   1e-9}  // 27 Air-Water Ud, delta - garbage
};

double Ca; // Ca = Mu*Ud/sigma
double Ca_mod; // Ca_mod = Mu*Umean/sigma
double Re; //Reynolds
double G;
double Umean, UMEAN;
double x_init = 2;
int maxlevel = 8;
int maxlevel_init = 13;
int minlevel = 5;
int LEVEL = 7;
int init_i = 1; // each restart it is equal to 1
int adapt_method = 1; // 0 - traditional, 1 - using limitation, 2 - using array for maxlevel
int snapshot_i = 100;
double fseps = 1e-3, ueps = 1e-2;
double TOLERANCE_P = 1e-5, TOLERANCE_V = 1e-5;
bool ellipse_shape = false, cylinder_shape = true;

int main (int argc, char * argv[]) {
    fprintf(ferr, "./a.out maxlevel bubcase adapt_method iter_fp lDomain dt_vtk snapshot_i\n");
//    maxruntime (&argc, argv);
    eta_s = 1e-5;
    TOLERANCE = 1e-6;
    NITERMIN = 1;
    NITERMAX = 100;
    DT = 1e-3;
    cylinder_shape = true;
//    relative_residual_poisson = true;
//    relative_residual_viscous = true;
    fs.refine = fs.prolongation = fraction_refine;
    f.refine = f.prolongation = fraction_refine;
    bubcase = 22;

    if (argc > 1)
        maxlevel = atoi (argv[1]);
    if (argc > 2)
        bubcase = atoi (argv[2]);

    if ((bubcase >= 8 && bubcase <= 12) || bubcase >= 18) { // Water
        Rho1 = 997, Rho2 = 1.204;
        Mu1 = 0.88e-3, Mu2 = 0.019e-3;
        Sigma = 72.8e-3;
        diam_tube = 514e-6;
        dt_vtk = 1e-1;
        lDomain = 20;
    }else { // Glycerol
        Rho1 = 1250, Rho2 = 1.204;
        Mu1 = 550e-3, Mu2 = 0.019e-3;
        Sigma = 63.4e-3;
        diam_tube = 494e-6;
        dt_vtk = 1e-2;
        lDomain = 10;
    }
    Ca = cases[bubcase].Ca; //Ca = Ud*Mu1/sigma
    Vd = cases[bubcase].Vd; // m^3
    UMEAN = cases[bubcase].Uc; // m/s

    if (argc > 3)
        adapt_method = atoi (argv[3]);
    if (argc > 4)
        iter_fp = atoi (argv[4]);
    if (argc > 5)
        lDomain = atof (argv[5]);
    if (argc > 6)
        dt_vtk = atof (argv[6]);
    if (argc > 7)
        snapshot_i = atoi (argv[7]);
    if (argc > 8)
        TOLERANCE_P = atof (argv[8]);
    if (argc > 9)
        TOLERANCE_V = atof (argv[9]);
    size (lDomain);
    origin (0., -L0/2., -L0/2.);
    init_grid (1 << LEVEL);

    deq = pow(6.0*Vd/pi, 1./3.);// 0.0005301091821 m
    dst = deq/diam_tube;// 1.0730955104
    rst = 0.5*dst;
    Vdst = (4./3.)*pi*cube(rst);
//    UMEAN = G*sq(0.5)/(8*Mu1);
    Ca_mod = Mu1*UMEAN/Sigma;
    Re = UMEAN*diam_tube*Rho1/Mu1;
    G = 32.0*Mu1*UMEAN/sq(diam_tube);

    if (ellipse_shape || dst < 0.9) {
        r_bub = min(rst, 0.4);
        l_bub = cube(rst) / sq(r_bub);
        ellipse_shape = true;
        cylinder_shape = false;
    }else if (cylinder_shape){
        r_bub = 0.45;
        l_bub = (Vdst - (4./3.)*pi*cube(r_bub))/(pi*sq(r_bub));
    } else{
        assert(false && "set shape");
    }
    x_init = 1.7*l_bub;

    fprintf(ferr,"BP:             eta_s=%g,     DT=%g\n"
                 "Solver:         NITERMIN=%d   NITERMAX=%d      TOLERANCE_P=%g TOLERANCE_V=%g TOLERANCE=%g  relative_residual_poisson=%d relative_residual_viscous=%d\n"
                 "OUTPUT:         dt_vtk=%g number of procs=%d\n"
                 "ADAPT:          minlevel=%d,  maxlevel=%d      adapt_meth=%d fseps=%g ueps=%g\n"
                 "Bubble case: %d\n"
                 "Properties(SI): Mu1=%g Mu2=%g Rho1=%g Rho2=%g  Sigma=%g G=%g UMEAN=%g\n"
                 "Apparatus:      diam_tube=%g  tube_length=%g\n"
                 "Bubble:         Vd=%g deq=%g  ellipse_shape=%d cylinder_shape=%d\n",
                 eta_s, DT,
                 NITERMIN, NITERMAX, TOLERANCE_P, TOLERANCE_V, TOLERANCE, relative_residual_poisson, relative_residual_viscous,
                 dt_vtk, npe(),
                 minlevel, maxlevel, adapt_method, fseps, ueps,
                 bubcase,
                 Mu1, Mu2, Rho1, Rho2, Sigma, G, UMEAN,
                 diam_tube, L0,
                 Vd, deq, ellipse_shape, cylinder_shape);
    // Dimensionless parameters:
    // Averaging on diam_tube=1 and Umean=1, Mu1=1 and Rho1=1 p' = p/(Rho1*UMEAN^2)
    G /= Rho1*sq(UMEAN)/diam_tube;
    Umean = 1;
    diam_tube /= diam_tube;
    RhoR = Rho1/Rho2;
    MuR = Mu1/Mu2;

    rho1 = 1.;// water
    rho2 = 1./RhoR; // air
    mu1 = 1./Re;
    mu2 = 1./(MuR*Re);
    f.sigma = 1./(Re*Ca_mod);
    fprintf(ferr,"Dimensionless Parameters: mu1=%g mu2=%g rho1=%g rho2=%g sigma=%g G=%g  Umean=%g\n"
                 "Dimensionless nums:       Re=%g  Ca=%g  Ca_mod=%g\n"
                 "Bubble:                   Vdst=%g dst=%g  rst=%g  r_bub=%g l_bub_cyl=%g l_bub=%g x_init=%g\n",
            mu1, mu2, rho1, rho2, f.sigma, G, Umean,
            Re, Ca, Ca_mod,
            Vdst, dst, rst, r_bub, l_bub, (ellipse_shape || dst < 0.9) ? l_bub : l_bub + 2*r_bub, x_init);
    int rank, psize, h_len;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    // get rank of this proces
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // get total process number
    MPI_Comm_size(MPI_COMM_WORLD, &psize);
    MPI_Get_processor_name(hostname, &h_len);
    printf("rank:%d size: %d at %s h_len %d\n", rank, psize, hostname, h_len);

    run();
}
//BCs

//Inflow
//u.n[left] = dirichlet((1 - fs[]));
u.n[left] = dirichlet(uexact(x,y,z)*(1 - fs[]));
p[left] = neumann(0);
pf[left] = neumann(0);
f[left] = dirichlet(1);
fs[left] = neumann(0);
//Outflow
u.n[right] = neumann(0);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);
f[right] = neumann(0);
fs[right] = neumann(0);
void geometry(scalar fs)
{
    //Channel
    vertex scalar phi[];
    foreach_vertex() {
        phi[] = sq(y) + sq(z) - sq(0.5);
    }
    boundary ({phi});
    fractions (phi, fs);
}

void bubble(scalar f)
{
    vertex scalar phi[];
    foreach_vertex() {
        if (ellipse_shape)
            phi[] = sq((x - x_init)/l_bub) + sq(y/r_bub) + sq(z/r_bub) - sq(1);
        else if (cylinder_shape)
        {
            phi[] = intersection( 0.5*l_bub - fabs(x - x_init),  sq(r_bub) - sq(y) - sq(z) );
            phi[] = union(phi[], sq(r_bub) - sq(fabs(x - x_init) - 0.5*l_bub) - sq(y) - sq(z));
            phi[] *= -1;
        }
    }
    boundary ({phi});
    fractions (phi, f);
}

int maXlevel(double x,double y, double z){
    double x0 = fabs(x - x_init - Umean*t);
    int n = ceil(max(0, 0.5*(x0/l_bub - 3)));
    return max(maxlevel-n, 10);
}

#define ADAPT_INIT {f, fs, u.x}
#define ADAPT_INIT_EPS {fseps, fseps, ueps}
#define ADAPT_INIT_MAXLEVEL {maxlevel, maxlevel, maxlevel-2}
event init (t = 0) {
    if (!restore (file = "restart")) {
        int it = 1;
        astats s;
        do {
            fprintf(ferr, "iteration=%d\n", it);
            count_cells(t, i);
            geometry(fs);
            bubble(f);
            foreach() {
                u.x[] = (1 - fs[])*uexact(x,y,z); //Init velocity
                u.y[] = 0;
                u.z[] = 0;
            }
            boundary((scalar *){u});
            if (adapt_method == 0)
                s = adapt_wavelet((scalar *) ADAPT_INIT, (double[]) ADAPT_INIT_EPS, maxlevel, minlevel);
            else if (adapt_method == 1)
                s = adapt_wavelet_limited((scalar *) ADAPT_INIT, (double []) ADAPT_INIT_EPS, maXlevel, minlevel);
            else if (adapt_method == 2)
                s = adapt_wavelet2((scalar *) ADAPT_INIT, (double[]) ADAPT_INIT_EPS, (int[]) ADAPT_INIT_MAXLEVEL, minlevel);
            fprintf(ferr, "Adaptation: nf=%d nc=%d\n", s.nf, s.nc);
            if (s.nf == 0  || it > 5) break;
            it++;
        } while(1);
//        event("vtk_file");
    }else{
        fprintf(ferr, "file is read with maxlevel_init=%d\n", grid->maxdepth);
    }
    maxlevel_init = grid->maxdepth;

}
event advection_term(i++){
    TOLERANCE = TOLERANCE_P;
}

event viscous_term(i++){
    TOLERANCE = TOLERANCE_V;
}

event projection(i++){
    TOLERANCE = TOLERANCE_P;
}

/**
## Counting bubbles

The number and sizes of bubbles is a useful statistics for atomisation problems.
This is not a quantity which is trivial to compute. The *tag()* function is
designed to solve this problem. Any connected region for which *f[] < 0.999*
(i.e. a bubble) will be identified by a unique "tag" value between 0 and *n-1*.
 */
event logfile (i +=1)
{
    double x_mean = 0, delta_min=1e+9, delta_mean=0, delta_max=-1e+9, r_min=1e+9;
    double x_min = 1e+9, x_max = -1e+9, length = 0, volume_clip = 0, volumeg = 0;
    double vel_bubx=0, vel_buby=0, vel_bubz=0;

    scalar m[];
    foreach() m[] = f[] < 0.999; // m is 0 and 1 array
    int n = tag (m); // m is modified filled be indices
    /**
    Once each cell is tagged with a unique droplet index, we can easily
    compute the volume *v* and position *b* of each droplet. Note that
    we use *foreach_leaf()* rather than *foreach()* to avoid doing a
    parallel traversal when using OpenMP. This is because we don't have
    reduction operations for the *v* and *b* arrays (yet).
     */
    double v[n];
    coord b[n];
    for (int j = 0; j < n; j++)
        v[j] = b[j].x = b[j].y = b[j].z = 0.;
    double vol_sum = 0;
    foreach_leaf()
        if (m[] > 0) {
            int j = m[] - 1; // j is index of a bubble
            v[j] += dv()*f[];
            vol_sum += v[j];
            coord p = {x,y,z};
            foreach_dimension()
                b[j].x += dv()*f[]*p.x;
        }
    /**
    When using MPI we need to perform a global reduction to get the
    volumes and positions of bubbles which span multiple processes. */

    #if _MPI
    MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    #endif

    /**
    Finally we output the volume and position of each bubble to
    standard output. */
    int i_vol_max; // index of the largest bubble
    double vol_max=-1e+9;
    for (int j = 0; j < n; j++){
        if (v[j] > vol_max) {
            vol_max = v[j];
            i_vol_max = j;
        }
        fprintf (stdout, "statistics: %d %g %d %g %g %g\n",
        i, t, j, v[j], b[j].x/v[j], b[j].y/v[j]);
    }
    fprintf(ferr, "i_vol_max= %d vol_max= %g\n", i_vol_max, vol_max);
    foreach(reduction(+:x_mean) reduction(+:volumeg)
            reduction(+:vel_bubx) reduction(+:vel_buby) reduction(+:vel_bubz)
            reduction(min:delta_min) reduction(min:x_min) reduction(max:x_max)
            ) {
        if (fs[]<1 && m[] == 1 + i_vol_max){ // inside of the largest bubble
            double dvtmp = (1.0 - f[])*dv(); // gas volume
            volumeg += dvtmp;
            x_mean   += x*dvtmp;// Along x
            vel_bubx += u.x[]*dvtmp; vel_buby += u.y[]*dvtmp; vel_bubz += u.z[]*dvtmp;//mean velocity of gas Ox,Oy,Oz
            if (f[] > 0 && f[] < 1) { // in bubble interface
                double r = sqrt(sq(y) + sq(z)) + Delta*(0.5 - f[]);
                if (0.5 - r < delta_min) {
                    delta_min = 0.5 - r;
                }
                if (x < x_min) x_min = x;
                if (x > x_max) x_max = x;
            }
        }
    }
    vel_bubx /= volumeg; vel_buby /= volumeg; vel_buby /= volumeg;
    x_mean /= volumeg; length = x_max - x_min;

    // Clip the largest bubble between xmin and x_mean
    foreach(reduction(+:volume_clip) reduction(max:delta_max) reduction(min:r_min)) {
        if ( f[] < 1 && m[] == 1 + i_vol_max){
            double r;
            if (x > x_min && x < x_mean && m[] == 1 + i_vol_max){ // look inside a clipped bubble
                volume_clip += (1.0 - f[])*dv(); // gas volume in the clipped volume
                if (f[] > 0 && 0.5 - r > delta_max) { // look at an interface of the clipped volume
                    r = sqrt(sq(y) + sq(z)) + Delta * (0.5 - f[]);
                    delta_max = 0.5 - r;
                }
            }
            if (f[] > 0 && x < x_mean) { // look at an left hemi-interface
                r = sqrt(sq(y) + sq(z));
                if (r < r_min) {
                    r_min = r;
//                    coord_tail.x = x;
//                    coord_tail.y = y;
//                    coord_tail.z = z;
                }
            }
        }
    }
    delta_mean = 0.5 - sqrt(volume_clip/(pi*(x_mean - x_min)));
        fprintf (ferr, "maxlevel= %d i= %d t= %g dt= %g volumeg= %g volume_clip= %g vel_bub= %g %g %g vel_bubx/U0-1= %g\n"
                   "x_min= %g x_mean= %g x_max= %g\n"
                   "delta_min= %g delta_mean(NOTE:x_mean_x_of_max_y)= %g delta_max= %g length= %g it_fp= %d\n",
            maxlevel, i, t, dt, volumeg, volume_clip, vel_bubx, vel_buby, vel_bubz, (vel_bubx/Umean - 1),
            x_min, x_mean, x_max,
            delta_min, delta_mean, delta_max, x_max - x_min, iter_fp);
    fprintf(ferr, "Ca\tUflow_m_s\tU_meanVT\tU_meanVT_m_s\tdelta_minVT\tdelta_meanVT\tdelta_maxVT\tmaxlevel\tlDomain\tdx\tN_per_delta\n");
    double dx_min = lDomain/pow(2., maxlevel);
   fprintf(ferr, "%8.5g\t%8.5g\t%8.5g\t%8.5g\t%8.5g\t%8.5g\t%8.5g\t%8d\t%8.5g\t%8.5g\t%8.5g\n",
            mu1*vel_bubx/f.sigma, UMEAN, vel_bubx, vel_bubx*UMEAN, delta_min, delta_mean, delta_max,
            maxlevel, lDomain, dx_min, delta_min/dx_min);
    if (i==0) fprintf(stdout, "t\tx_tail\tr_peak\tx_mean\tx_nose\tx_nose_ISC\tvolume\tUmeanV\t"
                              "delta_min\tdelta_mean\tdelta_max\tdelta_min_smooth\tdelta_max_smooth\n");
    fprintf (stdout, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
    t, x_min, 0.5 - delta_min, x_mean, x_max, 0.0, volumeg, vel_bubx,
    delta_min, delta_mean, delta_max, 0, 0);
    fflush(stdout);
}


event snapshot (i += snapshot_i)
//event snapshot (t += 1e-1)
{
    char name[80];
    scalar ppart[];
    scalar omega[];
    scalar l2[];
    sprintf(name, "dump-%04g",t);
    vorticity (u, omega);
    lambda2 (u, l2);
    foreach() ppart[] = pid();
    p.nodump = false;
    dump (file = name);
}

void exact(vector ue)
{
    foreach() {
        ue.x[] = (1 - fs[])*uexact(x,y,z);
        ue.y[] = 0;
        ue.z[] = 0;
    }
    boundary((scalar *){ue});
}
//event vtk_file (i += 1)
//event vtk_file (t += dt_vtk)
//{
//    char subname[80]; sprintf(subname, "tube_bp");
//    scalar l[], l2[], omega[];
//    foreach() l[] = level;
//    scalar np[]; foreach() np[] = pid();
//    vorticity (u, omega);
//    lambda2 (u, l2);
////    output_vtu_MPI( subname, (iter_fp) ? t + dt : 0, (scalar *) {p, fs, f, np, l}, (vector *) {u, a});
//    output_vtu_MPI( subname, (iter_fp) ? t + dt : 0, (scalar *) {p, fs, f, np, l, omega, l2}, (vector *) {u});
////    fprintf(ferr, "end:snapshot_vtk");
////    event("snapshot_vtk");
//}

int signnum(int x) {
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}

#define ADAPT_SCALARS {f, u.x, u.y, u.z}
#define ADAPT_EPS_SCALARS {fseps, ueps, ueps, ueps}
#define ADAPT_MAXLEVEL {maxlevel, max(maxlevel-2,10), max(maxlevel-2,10), max(maxlevel-2,10)}
//#define ADAPT_SCALARS {fs, f, u}
//#define ADAPT_EPS_SCALARS {fseps, fseps, ueps, ueps, ueps}
event adapt (i++)
{
    if (init_i % 100 == 0){
        maxlevel_init += signnum(maxlevel - maxlevel_init);
        fprintf(ferr, "Adaptation with maxlevel_init=%d", maxlevel_init);
    }
    double eps_arr[] = ADAPT_EPS_SCALARS;
    fprintf(ferr, "beginning adapt\n");
//    MinMaxValues (ADAPT_SCALARS, eps_arr);
    if (adapt_method == 0)
        adapt_wavelet ((scalar *) ADAPT_SCALARS, (double []) ADAPT_EPS_SCALARS, maxlevel = maxlevel_init, minlevel = minlevel);
    else if (adapt_method == 1)
//        adapt_wavelet_limited  ((scalar *) {f, u_mag}, (double []) {fseps, ueps}, maxlevel_init, minlevel);
        adapt_wavelet_limited  ((scalar *) ADAPT_SCALARS, (double []) ADAPT_EPS_SCALARS, maxlevel_init, minlevel);
    else if (adapt_method == 2)
        adapt_wavelet2((scalar *)ADAPT_SCALARS, (double []) ADAPT_EPS_SCALARS,(int []){maxlevel_init, maxlevel_init-1, maxlevel_init-2, maxlevel_init-2, maxlevel_init-2}, minlevel);

    fprintf(ferr, "ended adapt\n");
    count_cells(t, i);
    geometry(fs);
    double eps_arr2[] = {1, 1, 1, 1}; //u.x, u.y, u.z, p
    if (i % 10 == 0 || i < 5) MinMaxValues ({u, p}, eps_arr2);
}

event stop(t=L0/Umean);
