#define BRINKMAN_PENALIZATION 1
#define DEBUG_MINMAXVALUES
//#define DEBUG_BRINKMAN_PENALIZATION
#define DEBUG_MODE_POISSON
//#define DEBUG_OUTPUT_VTU_MPI
#define FILTERED
#define JACOBI 1

//#define PRINT_ALL_VALUES
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
//#define STOKES

scalar fs[];
scalar omega[];
scalar l2[];
scalar ppart[];
#include "grid/octree.h"
#include "../src_local/centered-weugene.h"
#include "two-phase.h"
#ifdef STOKES
    #include "navier-stokes/conserving.h"
#endif
#include "tension.h"
#include "../src_local/adapt_wavelet_limited.h"
#include "../src_local/adapt2.h"
#include "../src_local/utils-weugene.h"
#include "utils.h"
#include "lambda2.h"
#include "../src_local/output_vtu_foreach.h"
#include "maxruntime.h"

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
    double Ca;    // number of iterations
    double Vd;    // maximum residual before and after the iterations
    double Uc;    // sum of r.h.s.
    double Ud;    // number of relaxations
    double delta; // minimum level of the multigrid hierarchy
} Cases;
#define zv {0,0,0,0,0}
//         CA       V_d[l]    Uc[m/s]  Ud[m/s] delta*[-]
Cases cases[23]={
        zv, // 0 is empty
        zv,zv,zv,zv,zv,zv,zv, // 1-7 Air-glycerol
        {6.44e-4,  0.0349e-9, 0.0454, 0.0533, 0.1055},// 8  Air-Water delta?
        {8.42e-4,  0.0439e-9, 0.0593, 0.0697, 0.074}, // 9  Air-Water delta?
        {7.34e-4,  0.0602e-9, 0.0557, 0.0607, 0.027}, // 10 Air-Water delta?
        {6.64e-4,  0.0745e-9, 0.0534, 0.0550, 0.001}, // 11 Air-Water delta?
        {6.697e-4, 0.1893e-9, 0.0543, 0.0554, 0.006}, // 12 Air-Water
        zv,zv,zv,zv,zv, // 13-17 Air-glycerol
        {0.003, 0.1751e-9, 0.242, 0.261, 0.013},  // 18 Air-Water
        {0.008, 0.1715e-9, 0.666, 0.704, 0.023},  // 19 Air-Water
        {0.0098, 0.2208e-9, 0.757, 0.815, 0.025}, // 20 Air-Water
        {0.015, 0.1882e-9, 1.118 ,1.293, 0.039},  // 21 Air-Water
        {0.0230, 0.2179e-9, 1.580, 1.944, 0.054}  // 22 Air-Water
};

double Ca; // Ca = Mu*Ud/sigma
double Ca_mod; // Ca_mod = Mu*Umean/sigma
double Re; //Reynolds
double G;
double Umean;
double x_init = 2;
int maxlevel = 8;
int minlevel = 5;
int LEVEL = 6;
int adapt_method = 1; // 0 - traditional, 1 - using limitation, 2 - using array for maxlevel
int snapshot_i = 100;
double fseps = 1e-3, ueps = 1e-2;
double TOLERANCE_P = 1e-5, TOLERANCE_V = 1e-5;
bool ellipse_shape = false, cylinder_shape = true;
scalar un[];

int main (int argc, char * argv[]) {
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
    Umean = cases[bubcase].Uc; // m/s

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

    size(lDomain);
    origin(0., -L0/2., -L0/2.);
    init_grid(1 << LEVEL);

    deq = pow(6.0*Vd/pi, 1./3.);// 0.0005301091821 m
    dst = deq/diam_tube;// 1.0730955104
    rst = 0.5*dst;
    Vdst = (4./3.)*pi*cube(rst);
//    Umean = G*sq(0.5)/(8*Mu1);
    Ca_mod = Mu1*Umean/Sigma;
    Re = Umean*diam_tube*Rho1/Mu1;
    G = 32.0*Mu1*Umean/sq(diam_tube);

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
                 "Solver:         NITERMIN=%d   NITERMAX=%d      TOLERANCE=%g  relative_residual_poisson=%d relative_residual_viscous=%d\n"
                 "OUTPUT:         dt_vtk=%g number of procs=%d\n"
                 "ADAPT:          minlevel=%d,  maxlevel=%d      adapt_meth=%d fseps=%g ueps=%g\n"
                 "Bubble case: %d\n"
                 "Properties(SI): Mu1=%g Mu2=%g Rho1=%g Rho2=%g  Sigma=%g G=%g Umean=%g\n"
                 "Apparatus:      diam_tube=%g  tube_length=%g\n"
                 "Bubble:         Vd=%g deq=%g  ellipse_shape=%d cylinder_shape=%d resolve_capillary_effects=%d\n",
                 eta_s, DT,
                 NITERMIN, NITERMAX, TOLERANCE, relative_residual_poisson, relative_residual_viscous,
                 dt_vtk, npe(),
                 minlevel, maxlevel, adapt_method, fseps, ueps,
                 bubcase,
                 Mu1, Mu2, Rho1, Rho2, Sigma, G, Umean,
                 diam_tube, L0,
                 Vd, deq, ellipse_shape, cylinder_shape, resolve_capillary_effects);
    // Dimensionless parameters:
    // Averaging on diam_tube=1 and Umean=1, Mu1=1 and Rho1=1 p' = p/(Rho1*Umean^2)
    G /= Rho1*sq(Umean)/diam_tube;
    Umean /= Umean;
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
                 "Bubble:                   Vdst=%g dst=%g  rst=%g  r_bub=%g l_bub=%g x_init=%g\n",
            mu1, mu2, rho1, rho2, f.sigma, G, Umean,
            Re, Ca, Ca_mod,
            Vdst, dst, rst, r_bub, l_bub, x_init);
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
//            u.x[] = 1 - fs[]; //Init velocity
                u.x[] = (1 - fs[])*uexact(x,y,z); //Init velocity
                u.y[] = 0;
                u.z[] = 0;
                un[] = u.x[];
//            p[] = (1 - f[])*2.0*f.sigma/rst;
//            g.x[] = 0.5*((p[0] - p[-1])/Delta + (p[1] - p[0])/Delta);
            }
            boundary((scalar *){u, un});
            if (adapt_method == 0)
                s = adapt_wavelet((scalar *) ADAPT_INIT, (double[]) ADAPT_INIT_EPS, maxlevel, minlevel);
            else if (adapt_method == 1)
                s = adapt_wavelet_limited((scalar *) ADAPT_INIT, (double []) ADAPT_INIT_EPS, maXlevel, minlevel);
            else if (adapt_method == 2)
                s = adapt_wavelet2((scalar *) ADAPT_INIT, (double[]) ADAPT_INIT_EPS, (int[]) ADAPT_INIT_MAXLEVEL, minlevel);
            fprintf(ferr, "Adaptation: nf=%d nc=%d\n", s.nf, s.nc);
            if (s.nf == 0  || it > 8) break;
            it++;
        } while(1);
        event("vtk_file");
    }else{
        fprintf(ferr, "file is read\n");
    }
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

//
event logfile (i += 10) {
  double du = change (u.x, un);
  fprintf (ferr, "i= %d t= %g du= %g\n", i, t, du);
  if (i > 0 && du < 1e-6) //Convergence criteria
    return 1; //Stop the simulation
}

event logfile (i +=100)
{
    double avggas = L0*1 - normf(f).avg;
    scalar umag[];
    double xcg = 0, volume = 0, volumeg = 0, velamean = 0, velgx = 0, velgy = 0, dvtmp, min_r=0, mean_r=0, max_r=0, r=0;
    double length_min = 0, length_max = 0, length = 0;
    foreach(reduction(+:xcg) reduction(+:volume) reduction(+:volumeg)
            reduction(+:velgx) reduction(+:velgy) reduction(+:velamean)
            reduction(min:min_r) reduction(max:max_r) ) {
        if (fs[]<1){
            dvtmp = (1.0 - f[])*(1.0 - fs[])*dv(); // gas volume
            volumeg += dvtmp;//gas liquid
            volume += (1.0 - fs[])*dv();//channel volume
            umag[] = norm(u); // the length of u
            xcg   += x*dvtmp;// Along x
            velamean += umag[]*dvtmp;//mean velocity of gas
            velgx += u.x[]*dvtmp;//mean velocity of gas Ox
            velgy += u.y[]*dvtmp;//mean velocity of gas Oy
            if (f[] > 0 ) {
                r = sqrt(sq(y) + sq(z)) + Delta*(0.5 - f[]);
                if (r < min_r) min_r = r;
                if (r > max_r) max_r = r;
                length = x + Delta*(0.5 - f[]);
                if (length < length_min) length_min = length;
                if (length > length_max) length_max = length;
            }
        }
    }
    mean_r = 0.5*(min_r + max_r);
    xcg /= volumeg; velgx /= volumeg; velgy /= volumeg; velamean /= volumeg;
    //norm statu = normf_weugene(umag, fs); // outputs avg, rms, max, volume
    fprintf (ferr, "maxlevel= %d i= %d t= %g dt= %g avggas= %g velgx= %g valgy= %g velamean= %g velgx/U0-1= %g \n"
                   "xcg= %g thickness_min= %g thickness_mean= %g thickness_max= %g length= %g it_fp= %d\n",
            maxlevel, i, t, dt, avggas, velgx, velgy, velamean, (velgx/Umean - 1),
            xcg, 0.5 - max_r, 0.5 - mean_r, 0.5 - min_r, length_max - length_min, iter_fp);
}


//event movie (t += 0.1) {
//  view (fov = 22.4578, quat = {-0.707107,-0,-0,0.707107}, tx = -0.5, ty = 0.,
//  	bg = {0.3,0.4,0.6}, width = 600, height = 600, samples = 1);
//
//  clear();
//  squares("u.x", min=0, max=1.5, alpha = 0, n = {0,1,0});
//  save("movie_ux.mp4");
//}

event snapshot (i += snapshot_i)
//event snapshot (t += 1e-1)
{
    char name[80];
    sprintf(name, "dump-%04g",t);
    vorticity (u, omega);
    lambda2 (u, l2);
    foreach() ppart[] = pid();
    p.nodump = false;
    dump (file = name);
}

event snapshot_vtk (i += 10000)
{
    char name[80];
    sprintf(name, "dump-%04g", t);
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
//    fprintf(ferr, "l, ");
//    scalar l[]; foreach() l[] = level;
//    fprintf(ferr, "np, ");
//    scalar np[]; foreach() np[] = pid();
////    fprintf(ferr, "omega, ");
////    vorticity (u, omega);
////    fprintf(ferr, "lambda2, ");
////    lambda2 (u, l2);
////    fprintf(ferr, "output_vtu_MPI, ");
//
//    output_vtu_MPI( subname, (iter_fp) ? t + dt : 0, (scalar *) {p, fs, f, np, l}, (vector *) {u, a});
////    output_vtu_MPI( subname, (iter_fp) ? t + dt : 0, (scalar *) {p, fs, f, np, l, omega, l2}, (vector *) {u});
//    fprintf(ferr, "end:snapshot_vtk, ");
//    event("snapshot_vtk");
//}



#define ADAPT_SCALARS {f, u.x, u.y, u.z}
#define ADAPT_SCALARS_SMOOTHED {f, u.x, u.y, u.z}
#define ADAPT_EPS_SCALARS {fseps, ueps, ueps, ueps}
#define ADAPT_MAXLEVEL {maxlevel, max(maxlevel-2,10), max(maxlevel-2,10), max(maxlevel-2,10)}
//#define ADAPT_SCALARS {fs, f, u}
//#define ADAPT_EPS_SCALARS {fseps, fseps, ueps, ueps, ueps}
event adapt (i++)
{
    double eps_arr[] = ADAPT_EPS_SCALARS;
    fprintf(ferr, "beginning adapt\n");
    vorticity (u, omega);
//    MinMaxValues (ADAPT_SCALARS, eps_arr);
    if (adapt_method == 0)
        adapt_wavelet ((scalar *) ADAPT_SCALARS, (double []) ADAPT_EPS_SCALARS, maxlevel = maxlevel, minlevel = minlevel);
    else if (adapt_method == 1)
//        adapt_wavelet_limited  ((scalar *) {f, u_mag}, (double []) {fseps, ueps}, maXlevel, minlevel);
        adapt_wavelet_limited  ((scalar *) ADAPT_SCALARS, (double []) ADAPT_EPS_SCALARS, maXlevel, minlevel);
    else if (adapt_method == 2)
        adapt_wavelet2((scalar *)ADAPT_SCALARS, (double []) ADAPT_EPS_SCALARS,(int []){maxlevel, maxlevel-1, maxlevel-2, maxlevel-2, maxlevel-2}, minlevel);

    fprintf(ferr, "ended adapt\n");
    count_cells(t, i);
    geometry(fs);
    double eps_arr2[] = {1, 1, 1, 1}; //u.x, u.y, u.z, p
    if (i % 10 == 0 || i < 5) MinMaxValues ({u, p}, eps_arr2);
}

event stop(t=L0/Umean);
