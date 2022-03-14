#define DEBUG_MINMAXVALUES 1
#define EPS_MAXA 1
#define DEBUG_MULTIGRID
#define DEBUG_MODE_POISSON
#define PRINT_ALL_VALUES
#include "run.h"
#include "timestep.h"
#include "poisson-krylov.h"
//#include "poisson.h"
#include "utils.h"
#include "utils-weugene.h"
#include "output_vtu_foreach.h"

scalar u[], rhs[], Ap_result[], uexact[];
double w = 1.0;
double ueps = 1e-3;
double dtmax = 1e-21;
int minlevel = 2;
int maxlevel = 9;
int LEVEL = 7;

int main(int argc, char * argv[]) {
    DT = 1e+0;
    CFL = 0.5;
    TOLERANCE = 1e-8;
    size(2.0*pi);
    init_grid(1 << LEVEL);
    origin (-L0/2, -L0/2.);
//    periodic(right);
//    periodic(top);
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
//    for (maxlevel= 7; maxlevel <= 11; maxlevel++){
        run();
//    }
    return 0;
}

#define u_BC (sin(w*x*y))
#define du_BC (w*(y + x)*cos(w*x*y))// (w*y*cos + x)*cos(w*x*y))
#define d2u_BC (-sq(w)*(sq(y)+sq(x))*sin(w*x*y))
u[left] = dirichlet(u_BC);
u[right] = dirichlet(u_BC);
u[bottom] = dirichlet(u_BC);
u[top] = dirichlet(u_BC);

//u[left] = 0;
//u[right] = 0;
//u[bottom] = 0;
//u[top] = 0;

//u[right] = neumann (0);
//u[left]  = neumann (0);
//u[top]    = neumann (0);
//u[bottom] = neumann (0);

//rhs[left] = 0;
//rhs[right] = 0;
//rhs[bottom] = 0;
//rhs[top] = 0;

event init (t = 0) {
    foreach(){
        u[] = u_BC;
    }
    boundary((scalar *) {u});
    event("vtk_file");
}

event set_dtmax (i++,last) dtmax = DT;

event stability (i++, last) {
    dt = dtnext ( dtmax );
}

event step(i++){
//    diffusion (u, DT, r=rhs);
#if 1
    foreach(){
        rhs[] = 0;
//        rhs[] = d2u_BC;
    }
    boundary((scalar*) {rhs});
    poisson (u, rhs, unityf, zeroc, nrelax=4); //alpha=1, lambda=0
#else
    foreach(){
        rhs[] = 0;
    }
    boundary((scalar*) {rhs});

    struct Poisson p;
    p.a = u;
    p.b = rhs;
    p.alpha = unityf;
    p.lambda = zeroc;
    p.minlevel = 5;
    p.tolerance = 1e-7;
    p.nrelax = 100;

    if (!p.alpha.x.i)
        p.alpha = unityf;
    if (!p.lambda.i)
        p.lambda = zeroc;
    face vector alpha = p.alpha;
    scalar lambda = p.lambda;
    restriction ({alpha,lambda});


    residual ((scalar *){u}, (scalar *){rhs}, (scalar *){Ap_result}, &p); // res = b - Ay = -Ay
#endif
//    calcAy (u, rhs, p);

}

//#define ADAPT_SCALARS {u}
//#define ADAPT_EPS_SCALARS {ueps}

//event adapt (i++){
//    double eps_arr[] = ADAPT_EPS_SCALARS;
//    MinMaxValues((scalar *) ADAPT_SCALARS, eps_arr);
//    adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
//
//}

event vtk_file (i+=1)
{
    char name[300];
    sprintf (name, "vtk_krylov");
    scalar l[];
    foreach() {
        l[] = level;
        uexact[] = u_BC;
    }
    boundary((scalar *){l, uexact});

    output_vtu_MPI(name, (iter_fp) ? t + dt : 0, list = (scalar *) {u, uexact, rhs, Ap_result, l, residual_of_p});

}





event stop_end(i = 10){
    scalar du[];
    foreach(){
        du[] = fabs(u[] - uexact[]);
    }
    boundary((scalar *){du});
    double eps_arr[] = {1};
    MinMaxValues((scalar *){du}, eps_arr);
    //count_cells(t, i);
    fprintf(fout, "L2 h=%g L2Err=%g Nc=%ld\n", L0/pow(2.0, maxlevel), eps_arr[0], perf.tnc);

}
#if TRACE > 1
event profiling (i++) {
  static FILE * fp = fopen ("profiling", "w");
  trace_print (fp, 1);
}
#endif
