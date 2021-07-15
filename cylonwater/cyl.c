#define BRINKMAN_PENALIZATION 1
#define DEBUG_MINMAXVALUES
#define DEBUG_BRINKMAN_PENALIZATION
#define DEBUG_MODE_TENSION
#define DEBUG_MODE_POISSON
//#define DEBUG_OUTPUT_VTU_MPI
#define FILTERED
#define JACOBI 1
#define DAMP_CAPILLARY_WAVE
#define CURVATURE_CORR
#define CORRECT_UF_FLUXES 1
#define RELATIVE_RES
#define EPS_MAXA 2


#ifdef DEBUG_MODE_TENSION
    scalar f_hat[];
#endif
scalar which_meth[];
face vector fs_face[];
//#include "grid/octree.h"
#include "curvature_partstr.h"
#include "centered-weugene.h"
#include "three-phase-weugene.h"
#include "tension.h"
#include "adapt_wavelet_limited.h"
#include "adapt2.h"
#include "utils-weugene.h"
#include "rvachev.h"
#include "output_vtu_foreach.h"
#if dimension > 2
    #include "lambda2.h"
#endif
int maxlevel = 8;
int minlevel = 6;
double U0=1, rhol=1, sig=0.0005, Lchar=1, mul=1;
double Rrho=0.1, Rmu=0.1, rad=0.0625, maxDT, ys0=0.0009765625;
double RE, CA;
double rhoeps=1e-8, feps=1e-8, ueps=1e-3;
double TOLERANCE_P = 1e-8, TOLERANCE_V = 1e-7;
char subname[80];
#if MOVE==1
    static coord vel_s = {-1, 0, 0};
    double xs0 = 0.3125;
#else
    static coord vel_s = {0, 0, 0};
    double xs0 = -0.3125;
#endif

/**
The domain is the periodic unit square centered on the origin. */

int main(int argc, char * argv[])
{
    fprintf(ferr, "./a.out RE CA rad Rrho Rmu maxlevel for %dD problem\n", dimension);
    RE = U0*2*rad*rhol/mul;
    CA = U0*mul/sig;
    TOLERANCE = 1e-6;
    maxDT = 1e-4;
    m_bp = 0.1; // Brinkman layer resolution
    NITERMIN = 5;
    NITERMAX = 100;
    DT = 1e-6;
    CFL = 0.4;
    ignore_tension_CFL = false;
    CFL_SIGMA = 0.5;// it is ignored if ignore_tension_CFL=true

    if (argc > 1) {
        RE = atof(argv[1]);
    }
    if (argc > 2) {
        CA = atof(argv[2]);
    }
    if (argc > 3) {
        rad = atof(argv[3]);
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
    if (argc > 7) {
        strcpy(subname, argv[7]);
    }

    size (1.0);
    origin (-0.5*L0, -0.5*L0, -0.5*L0);
    init_grid( 1 << minlevel );
    periodic(right);

    rho1 = 1.; rho2 = rho1*Rrho; rho3 = 1.5*max(rho1, rho2);
    mu1 = 1./RE; mu2 = mu1*Rmu; mu3 = mu1*rho3/rho1;
    f.sigma = 1./(RE*CA);
#if HYDROPHOBIC == 1
    ys0 = -L0/pow(2, maxlevel+2);
#else
    ys0 = L0/pow(2, maxlevel+2);
#endif
    fprintf(ferr, "maxlevel=%d tol=%g NITERMAX=%d NITERMAX=%d\n"
                  "RE=%g CA=%g rad=%g rho1/rho2=%g mu1/mu2=%g\n"
                  "mu1=%g mu2=%g rho1=%g rho2=%g sigma=%g\n"
                  "xs0=%g, ys0=%g\n",
            maxlevel, TOLERANCE, NITERMIN, NITERMAX,
            RE, CA, rad, Rrho, Rmu,
            mu1, mu2, rho1, rho2, f.sigma,
            xs0, ys0);

    const vector U_sol[] = {vel_s.x, vel_s.y, vel_s.z};
    target_U = U_sol;

    fs.refine = fs.prolongation = fraction_refine;
    f.refine = f.prolongation = fraction_refine;
    boundary((scalar *){f, fs});
#ifdef _MPI
    int rank, psize, h_len;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    // get rank of this proces
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // get total process number
    MPI_Comm_size(MPI_COMM_WORLD, &psize);
    MPI_Get_processor_name(hostname, &h_len);
    printf("rank:%d size: %d at %s h_len %d\n", rank, psize, hostname, h_len);
#endif
    DumpCsvInit();
    partstr_conf.nohf = false;
    run();
    DumpCsvFin();
}

//Inflow

#if dimension > 2
    u.r[bottom] = neumann(0);
    u.r[top] = neumann(0);
    u.r[back] = neumann(0);
    u.r[front] = neumann(0);

    u.n[back] = dirichlet(0);
    u.t[back] = neumann(0);
    p[back] = neumann(0);
    pf[back] = neumann(0);
    f[back] = neumann(0);
    fs[back] = neumann(0);

    u.n[front] = dirichlet(0);
    u.t[front] = neumann(0);
    p[front] = neumann(0);
    pf[front] = neumann(0);
    f[front] = neumann(0);
    fs[front] = neumann(0);
#endif


u.n[bottom] = dirichlet(0);
u.t[bottom] = neumann(0);
p[bottom] = neumann(0);
pf[bottom] = neumann(0);
f[bottom] = neumann(0);
fs[bottom] = neumann(0);

u.n[top] = dirichlet(0);
u.t[top] = neumann(0);
p[top] = neumann(0);
pf[top] = neumann(0);
f[top] = neumann(0);
fs[top] = neumann(0);

scalar divu[];
#define tmpposx (xs0 + vloc*t)
#define outOfBox ((positionX < X0) || (positionX > X0 +L0))
double posx (double t, double vloc){
    double positionX = tmpposx;
    double signvc = (vloc > 0) ? 1 : (vloc < 0)? -1 : 0;
    while(outOfBox){
        positionX -= signvc*L0;
    }
    return positionX;
}
void soild_fs(scalar fs, face vector fs_face, double t){
    vertex scalar phi[];
    foreach_vertex() {
        phi[] = -HUGE;
        for (int xi=-L0; xi <=L0; xi +=L0) {
            phi[] = max(phi[], (sq(rad) - sq(x - posx(t, vel_s.x) - xi) - sq(y) - sq(z)));
        }
    }
    boundary ({phi});
    fractions (phi, fs, fs_face);
}

#define f1(x,y,z) (-(y) + ys0)
#define f2(x,y,z) (sq(rad) - sq((x) - posx(t, vel_s.x) - xi) - sq(y) - sq(z))
#define f3(x,y,z) (sq(0.2*rad) - sq(sqrt(sq((x) - posx(t, vel_s.x) - xi) + sq(y)) - rad) - sq(z))


void bubbles (scalar f)
{
    vertex scalar phi[];
    face vector ff[];
    double a0=1, a1=1, a2=1;
    foreach_vertex() {
        phi[] = f1(x,y,z);
        for (int xi=-L0; xi <=L0; xi +=L0) {
#if HYDROPHOBIC == 1
            phi[] = CBRsubstraction(phi[], f2(x - xi,y,z), f3(x - xi,y,z), -0.001, 1, 1, 0.05); // yama
//            phi[] = min(phi[], -f2(x - xi,y,z)); // yama
#else
            phi[] = CBRunion(phi[], f2(x - xi,y,z), f3(x - xi,y,z), 0.001, 1, 1, 0.05); // bugor
//            phi[] = max(phi[],  f2(x - xi,y,z)); // bugor
#endif
        }
    }
    boundary ({phi});
    fractions (phi, f, ff);
}


event init (t = 0) {
    if (!restore (file = "restart")) {
        int it = 0;
        do {
            soild_fs (fs, fs_face, 0);
            bubbles(f);
        }while (adapt_wavelet({f, fs}, (double []){feps, feps}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
        refine( (fabs(y) < L0/pow(2, minlevel)) && level < maxlevel );
        refine( ((sq(x - xs0) + sq(y) + sq(z) < sq(1.1 * rad)) && (sq(x - xs0) + sq(y) + sq(z) > sq(0.9 * rad))) && level < maxlevel );
        soild_fs (fs, fs_face, 0);
        bubbles(f);

        foreach() {
#if MOVE==1
            foreach_dimension() u.x[] = fs[]*target_U.x[];
#else
            u.x[] = (1 - fs[])*U0;
#endif
        }
        boundary((scalar *){u});
        event("vtk_file");
    }
}

event set_dtmax (i++) {
    RELATIVE_RES_TOLERANCE = 0.04; //fabs(res1 - res2)/(res1 + 1e-30) < RELATIVE_RES_TOLERANCE
    DT *= 1.05;
    DT = min(DT, maxDT);
    fprintf(ferr, "TIMEMAX: tnext= %g, t=%g, DT=%g, dt=%g\n", tnext, t, DT, dt);
}

event advection_term (i++){
    TOLERANCE = TOLERANCE_P;
//    soild_fs (fs, fs_face, t + dt/2.);
}

event viscous_term (i++){
    TOLERANCE = TOLERANCE_V;
//    soild_fs (fs, fs_face, t + dt);
}

event projection(i++){
    TOLERANCE = TOLERANCE_P;
#if MOVE==1
    for (int xi=-L0; xi <=L0; xi +=L0) {
        if (fabs(fabs(posx(t, vel_s.x) + xi) - fabs(xs0)) < 2*U0*dt) {
            event("vtk_file");
            break;
        }
    }
#endif
}

event lof_file(i += 20){
    double V1=0;
    double u_mean_x=0, u_max_x=0, u1_mean_x=0, u1_max_x=0;
    foreach(reduction(+:V1) reduction(+:u_mean_x) reduction(max:u_max_x)
            reduction(+:u1_mean_x) reduction(max:u1_max_x)
    ){
        V1 += f[]*dv();
        u_mean_x += u.x[]*dv();
        u1_mean_x += f[]*u.x[]*dv();
        if (fabs(u.x[]) > u_max_x)
            u1_max_x = u.x[];
        if (fabs(u.x[]*f[]) > u1_max_x)
            u1_max_x = u.x[]*f[];
    }
    fprintf(ferr, "logs: V1: %g u_mean_x: %g u1_mean_x: %g  u_max_x: %g u1_max_x: %g\n",
                         V1, u_mean_x, u1_mean_x, u_max_x, u1_max_x);
}

//Output
//event vtk_file (i += 10){
event vtk_file (t += 0.01){
    scalar l[], omega[], divu[], l2[];
    face vector uf_corr[];
    foreach() {
        l[] = level;
        divu[] = 0;
        foreach_dimension() divu[] += (uf.x[1]-uf.x[])/Delta;
    }
    boundary((scalar *) {l, divu});
    foreach_face() {
        uf_corr.x[] = fs_face.x[]*(uf.x[] - face_value(target_U.x, 0));
    }
    boundary((scalar *){uf_corr});
    brinkman_correction_uf (uf);
    vorticity (u, omega);
    #if dimension > 2
        lambda2 (u, l2);
    #endif
    output_vtu_MPI(subname, (iter_fp) ? t + dt : 0,
               list = (scalar *) {p, divu, fs, f, f_hat, l, rhov, omega, l2, which_meth, my_kappa},
               vlist = (vector *) {u, g, uf, uf_corr, a, dbp, total_rhs, divtauu, fs_face});
}

#if dimension > 2
    #define ADAPT_SCALARS {rhov, fs, u.x, u.y, u.z}
    #define ADAPT_EPS_SCALARS {rhoeps, feps, ueps, ueps, ueps}
#else
    #define ADAPT_SCALARS {rhov, fs, u.x, u.y}
    #define ADAPT_EPS_SCALARS {feps, feps, ueps, ueps}
#endif
event adapt (i++){
    double eps_arr[] = ADAPT_EPS_SCALARS;
    MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
    soild_fs (fs, fs_face, t + dt); // be careful with dt !!!
}

event snapshot (t += 0.1) {
      char name[80];
      sprintf (name, "dump_%s_%dD-%g", subname, dimension, t);
      dump (name);
}

event stop(t = 25*fabs(xs0)/U0);