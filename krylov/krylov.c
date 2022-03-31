#define DEBUG_MINMAXVALUES 1
#define EPS_MAXA 1
#define DEBUG_MULTIGRID
#define DEBUG_MODE_POISSON
#define PRINT_ALL_VALUES
const face vector unityfn[] = {-1.,-1.,-1.};
scalar z0[], z1[], res0[], res1[], p0[], Ap[];
scalar da[];
#include "run.h"
#include "timestep.h"
#include "poisson-krylov.h"
//#include "poisson.h"
#include "utils.h"
#include "utils-weugene.h"
#include "output_htg.h"
//#include "output_vtu_foreach.h"
scalar u[], rhs[], Ap_result[], uexact[];
double w = 1.0;
double ueps = 1e-3;
double dtmax = 1e-21;
int minlevel = 2;
int maxlevel = 9;
int LEVEL = 6;

int main(int argc, char * argv[]) {
    DT = 1e+0;
    CFL = 0.5;
    TOLERANCE = 1e-12;
    size(2.0);
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

//#define u_BC (sin(2*pi*x)*sin(3*pi*y))
//#define d2u_BC (-13*sq(pi)*u_BC)
#define u_BC (sin(2*pi*x)*sin(3*pi*y))
#define d2u_BC (13*sq(pi)*u_BC)
//#define u_BC (sin(x)*sin(y)*(sq(x) - sq(pi)))
//#define du_BC (w*(y + x)*cos(w*x*y))// (w*y*cos + x)*cos(w*x*y))
//#define d2u_BC (-sq(w)*(sq(y)+sq(x))*sin(w*x*y))

u[left] = dirichlet(u_BC);
u[right] = dirichlet(u_BC);
u[bottom] = dirichlet(u_BC);
u[top] = dirichlet(u_BC);

event init (t = 0) {
    foreach(){
        u[] = 0.5*noise();
//        u[] = u_BC;
    }
    event("report");
}

event set_dtmax (i++,last) dtmax = DT;

event stability (i++, last) {
    dt = dtnext ( dtmax );
}

event step(i++){
    foreach(){
        rhs[] = d2u_BC;
    }
    poisson (u, rhs, unityfn, zeroc); //alpha=1, lambda=0: div( unityfn*grad u ) = rhs

}

//#define ADAPT_SCALARS {u}
//#define ADAPT_EPS_SCALARS {ueps}

//event adapt (i++){
//    double eps_arr[] = ADAPT_EPS_SCALARS;
//    MinMaxValues((scalar *) ADAPT_SCALARS, eps_arr);
//    adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
//
//}

//event vtk_file (i+=1)
//{
//    char name[300];
//    sprintf (name, "vtk_krylov");
//    scalar l[];
//    foreach() {
//        l[] = level;
//        uexact[] = u_BC;
//    }
//
//    output_vtu_MPI(name, (iter_fp) ? t + dt : 0, list = (scalar *) {u, uexact, rhs, Ap_result, level, residual_of_p, res0,res1,p0,Ap}, fvlist = (vector*){g});
//
//}


event report(i++){
    char path[] = "res"; // no slash at the end!!
    char prefix[] = "data";
    scalar l[], laplaceU[];
    foreach() {
        l[] = level;
        uexact[] = u_BC;
//        laplaceU[] = (u[-1] - 2 * u[] + u[1]) / sq(Delta) + (u[0,-1] - 2 * u[] + u[0,1]) / sq(Delta);
        laplaceU[] = 0;
        foreach_dimension()
            laplaceU[] += (u[-1] - 2 * u[] + u[1]) / sq(Delta);
    }
    output_htg(path, prefix,  (iter_fp) ? t + dt : 0, (scalar *) {u, uexact, rhs, Ap_result, l, res0, res1, p0, Ap, laplaceU});
}
#if TRACE > 1
    event profiling (i += 20) {
      static FILE * fp = fopen ("profiling", "w");
      trace_print (fp, 1);
    }
#endif

event stop_end(i = 0){
    scalar du[];
    foreach(){
        du[] = fabs(u[] - uexact[]);
    }
    double eps_arr[] = {1};
    MinMaxValues((scalar *){du}, eps_arr);
    fprintf(fout, "L2 h=%g L2Err=%g Nc=%ld\n", L0/pow(2.0, maxlevel), eps_arr[0], perf.tnc);

}