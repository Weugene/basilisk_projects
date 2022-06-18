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
scalar omega[];
scalar l2[];
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
#include "output_vtu_foreach.h"
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
int minlevel = 5;
int LEVEL = 7;
int adapt_method = 1; // 0 - traditional, 1 - using limitation, 2 - using array for maxlevel
int snapshot_i = 100;
double fseps = 1e-3, ueps = 1e-2;
double TOLERANCE_P = 1e-5, TOLERANCE_V = 1e-5;
char dumpname[150];
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
    lDomain = 30;
    sprintf (dumpname, "restart");

    if (argc > 1)
        maxlevel = atoi (argv[1]);
    if (argc > 2)
        bubcase = atoi (argv[2]);
    if (argc > 3)
        lDomain = atof (argv[3]);
    if (argc > 4)
        sprintf (dumpname, "%s", argv[4]);


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

    fprintf(ferr,"BP:             eta_s=%g,     DT=%g, dumpname=%s\n"
                 "Solver:         NITERMIN=%d   NITERMAX=%d      TOLERANCE=%g  relative_residual_poisson=%d relative_residual_viscous=%d\n"
                 "OUTPUT:         dt_vtk=%g number of procs=%d\n"
                 "ADAPT:          minlevel=%d,  maxlevel=%d      adapt_meth=%d fseps=%g ueps=%g\n"
                 "Bubble case: %d\n"
                 "Properties(SI): Mu1=%g Mu2=%g Rho1=%g Rho2=%g  Sigma=%g G=%g UMEAN=%g\n"
                 "Apparatus:      diam_tube=%g  tube_length=%g\n"
                 "Bubble:         Vd=%g deq=%g  ellipse_shape=%d cylinder_shape=%d\n",
            eta_s, DT, dumpname,
            NITERMIN, NITERMAX, TOLERANCE, relative_residual_poisson, relative_residual_viscous,
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

event init (t = 0) {
    fprintf(ferr, "L before: %g\n", L0);
    if (!restore (file = "restart")) {
        fprintf(ferr, "file can\'t be read\n");
    }else{
        fprintf(ferr, "file is read\n");
        fprintf(ferr, "L after: %g\n", L0);
    }
}

#define ADAPT_SCALARS {u.x}
#define ADAPT_EPS_SCALARS {ueps}

event adapt_custom (i += 1)
{
    double eps_arr[] = ADAPT_EPS_SCALARS;
    fprintf(ferr, "beginning adapt maxlevel=%d\n", maxlevel);
    vorticity (u, omega);
    fprintf(ferr, "vorticity done\n");
    MinMaxValues((scalar *) ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);

    fprintf(ferr, "ended adapt\n");
    count_cells(t, i);
    geometry(fs);
}

event snapshot (i += 1)
{
    char name[80];
    sprintf(name, "restart_maxlevel-%d-%04g", maxlevel, t);
    fprintf(ferr, "snapshot: %s", name);
    vorticity (u, omega);
    lambda2 (u, l2);
    p.nodump = false;
    dump (file = name);
    exit(1);
}

event stop(t=L0/Umean);
