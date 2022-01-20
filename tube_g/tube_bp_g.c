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
face vector fs_face=zerof;
scalar omega[];
scalar l2[];
#include "grid/octree.h"
#include "centered-weugene.h"
#include "two-phase.h"
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
double Vdst, deq_st = 0.2, req_st = 0.1, r_bub, l_bub;
double RhoR, MuR;
double Ggrav;
double dt_vtk = 0.1;
double lDomain;
int bubcase;

double Ca; // Ca = Mu*Ud/sigma
double Re; //Reynolds number
double Bo; //Buoyancy number
double G;
double x_init = 2;
int maxlevel = 8;
int minlevel = 5;
int LEVEL = 7;
int adapt_method = 0; // 0 - traditional, 1 - using limitation, 2 - using array for maxlevel
int snapshot_i = 100;
double fseps = 1e-3, ueps = 1e-2;
double TOLERANCE_P = 1e-7, TOLERANCE_V = 1e-7;
bool ellipse_shape = false, cylinder_shape = true;
char subname[150];

int main (int argc, char * argv[]) {
    fprintf(ferr, "./a.out maxlevel bubcase adapt_method iter_fp lDomain dt_vtk snapshot_i\n");
//    maxruntime (&argc, argv);
    eta_s = 1e-5;
    TOLERANCE = 1e-7;
    NITERMIN = 1;
    NITERMAX = 100;
    DT = 1e-3;
    cylinder_shape = true;
//    relative_residual_poisson = true;
//    relative_residual_viscous = true;
    fs.refine = fs.prolongation = fraction_refine;
    f.refine = f.prolongation = fraction_refine;


    if (argc > 1)
        maxlevel = atoi(argv[1]);
    if (argc > 2)
        Re = atof(argv[2]);
    if (argc > 3)
        Ca = atof(argv[3]);
    if (argc > 4)
        Bo = atof(argv[4]);
    if (argc > 5)
        RhoR = atof(argv[5]);
    if (argc > 6)
        MuR = atof(argv[6]);
    if (argc > 7)
        deq_st = atof(argv[7]);
    if (argc > 8)
        lDomain = atof(argv[8]);
    if (argc > 9)
        adapt_method = atoi(argv[9]);
    if (argc > 10)
        iter_fp = atoi(argv[10]);
    if (argc > 11)
        dt_vtk = atof(argv[11]);
    if (argc > 12)
        snapshot_i = atoi(argv[12]);
    if (argc > 13)
        strcpy(subname, argv[13]);

    rho1 = 1.;// water
    rho2 = 1./RhoR; // air
    mu1 = 1./Re;
    mu2 = 1./(MuR*Re);
    f.sigma = 1./(Re*Ca);
    Ggrav = -4.*Bo/(Re*Ca);
    size (lDomain);
    origin (0., -L0/2., -L0/2.);
    init_grid (1 << LEVEL);
    periodic(top);
    periodic(front);
    Vdst = pi*cube(deq_st)/6.0;
    req_st = 0.5*deq_st;


    if (ellipse_shape || deq_st < 0.9) {
        r_bub = min(req_st, 0.4);
        l_bub = cube(req_st) / sq(r_bub);
        ellipse_shape = true;
        cylinder_shape = false;
    }else if (cylinder_shape){
        r_bub = 0.45;
        l_bub = (Vdst - (4./3.)*pi*cube(r_bub))/(pi*sq(r_bub));
    } else{
        assert(false && "set shape");
    }
    x_init = max(0.5, 1.7*l_bub);
    fprintf(ferr,"BP:             eta_s=%g,     DT=%g\n"
                 "Solver:         NITERMIN=%d   NITERMAX=%d      TOLERANCE=%g  relative_residual_poisson=%d relative_residual_viscous=%d\n"
                 "OUTPUT:         dt_vtk=%g number of procs=%d snapshot_i=%d subname=%s\n"
                 "ADAPT:          minlevel=%d,  maxlevel=%d      adapt_meth=%d fseps=%g ueps=%g\n"
                 "Geometry:       tube_length=%g, Vdst=%g deq_st=%g, r_bub=%g, l_bub=%g, ellipse_shape=%d, cylinder_shape=%d, x_init=%g\n"
                 "Dimensionless nums:       Re=rho1*diam*Umean/mu1=%g,  Ca=mu1*Umean/sigma=%g,  Bo=rho1*g*R^2/sigma=%g, RhoR=%g, MuR=%g, lDomain=%g\n"
                 "Dimensionless Parameters: mu1=%g, mu2=%g, rho1=%g, rho2=%g, sigma=%g, Ggrav=%g\n",
            eta_s, DT,
            NITERMIN, NITERMAX, TOLERANCE, relative_residual_poisson, relative_residual_viscous,
            dt_vtk, npe(), snapshot_i, subname,
            minlevel, maxlevel, adapt_method, fseps, ueps,
            L0, Vdst, deq_st, r_bub, l_bub, ellipse_shape, cylinder_shape, x_init,
            Re, Ca, Bo, RhoR, MuR, lDomain,
            mu1, mu2, rho1, rho2, f.sigma, Ggrav);

#if _MPI
    int rank, psize, h_len;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    // get rank of this proces
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // get total process number
    MPI_Comm_size(MPI_COMM_WORLD, &psize);
    MPI_Get_processor_name(hostname, &h_len);
    printf("rank:%d size: %d at %s h_len %d\n", rank, psize, hostname, h_len);
#endif
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
    double x0 = fabs(x - x_init - t);
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

event acceleration (i++) {
    if (Ggrav){
        face vector av = a;
        foreach_face(y) av.y[] += Ggrav;
    }
}

event projection(i++){
    TOLERANCE = TOLERANCE_P;
}

event snapshot (i += snapshot_i)
//event snapshot (t += 1e-1)
{
    if (i>0){
        char name[80];
        scalar ppart[];
        sprintf(name, "dump-%04g",t);
        vorticity (u, omega);
        lambda2 (u, l2);
        foreach() ppart[] = pid();
        p.nodump = false;
        dump (file = name);
    }
}

event snapshot_vtk (i += 100000)
{
    char name[300];
    sprintf (name, "dump_%s-%04g", subname, t);
    dump (name);
}


//event vtk_file (i += 100)
event vtk_file (t += dt_vtk)
{
    char name[300];
    sprintf (name, "vtk_%s", subname);
    scalar l[]; foreach() l[] = level;
    vorticity (u, omega);
    lambda2 (u, l2);
    output_vtu_MPI( name, (iter_fp) ? t + dt : 0, (scalar *) {p, fs, f, l, rhov, omega, l2}, (vector *) {u, a, g});
}



#define ADAPT_SCALARS {f, u.x, u.y, u.z}
#define ADAPT_EPS_SCALARS {fseps, ueps, ueps, ueps}
#define ADAPT_MAXLEVEL {maxlevel, max(maxlevel-2,10), max(maxlevel-2,10), max(maxlevel-2,10)}
//#define ADAPT_SCALARS {fs, f, u}
//#define ADAPT_EPS_SCALARS {fseps, fseps, ueps, ueps, ueps}
event adapt (i++)
{
    double eps_arr[] = ADAPT_EPS_SCALARS;
    fprintf(ferr, "beginning adapt\n");
    vorticity (u, omega);
    MinMaxValues (ADAPT_SCALARS, eps_arr);
    if (adapt_method == 0)
        adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
    else if (adapt_method == 1)
//        adapt_wavelet_limited  ((scalar *) {f, u_mag}, (double []) {fseps, ueps}, maXlevel, minlevel);
        adapt_wavelet_limited  ((scalar *) ADAPT_SCALARS, eps_arr, maXlevel, minlevel);
    else if (adapt_method == 2)
        adapt_wavelet2((scalar *)ADAPT_SCALARS, eps_arr,(int []){maxlevel, maxlevel-1, maxlevel-2, maxlevel-2, maxlevel-2}, minlevel);

    fprintf(ferr, "ended adapt\n");
    count_cells(t, i);
    geometry(fs);
    // Statistics
    double eps_arr2[] = {1, 1, 1, 1}; //u.x, u.y, u.z, p
    if (i % 10 == 0 || i < 5) MinMaxValues ({u, p}, eps_arr2);
}

event stop(t=L0);
