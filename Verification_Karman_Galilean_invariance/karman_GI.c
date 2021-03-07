#define BRINKMAN_PENALIZATION 1
#define DEBUG_BRINKMAN_PENALIZATION 1
//#define DEBUG_OUTPUT_VTU_MPI
#define DEBUG_MINMAXVALUES
#define DEBUG_MODE
//#define MOVING_ON
#define FILTERED
#define JACOBI 1
#define EPS_MAXA 2
#define RELATIVE_RESIDUAL
#define MODIFIED_CHORIN 0
#ifdef DEBUG_MODE
    scalar divutmpAfter[];
    scalar mod_du_dx[];
    vector conv_term[];
#endif
#include "../src_local/centered-weugene.h"
//#include "../src_local/three-phase-weugene.h"//rho3, mu3
#include "view.h"
#include "../src_local/output_vtu_foreach.h"
#include "tracer.h"
scalar f[];
scalar * tracers = {f};

int maxlevel;
int minlevel;
double eps = 1e-2;
double Ldomain = 1.0, RE=40.;
double xx0 = 0, rad = 0.1;
scalar fs[], omega[];
double signvc = -1;
#ifdef MOVING_ON
    coord vc = {-1.0, 0.0, 0.0};
    double vin = 0;
#else
    coord vc = {0.0, 0.0, 0.0};
    double vin = 1;
#endif
/**
The domain is the periodic unit square centered on the origin. */

/**
The fluid is injected on the left boundary with a unit velocity. The
tracer is injected in the lower-half of the left boundary. An outflow
condition is used on the right boundary. */
#ifdef MOVING_ON
u.n[left]  = neumann(0.);
u.t[left]  = neumann(0.);
uf.n[left] = neumann(0.);
uf.t[left] = neumann(0.);
#else
u.n[left]  = dirichlet(vin);
u.t[left]  = dirichlet(0.);
uf.n[left] = vin;
uf.t[left] = 0;
#endif
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
//f[left]   = dirichlet(y < 0);

u.n[right] = neumann(0);
u.t[right] = neumann(0);
uf.n[right]  = neumann(0.);
uf.t[right]  = neumann(0.);
#ifdef MOVING_ON
p[right]   = neumann(0.);
pf[right]  = neumann(0.);
#else
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);
#endif
#define tmpposx (xx0 + vc.x*t)
#define outOfBox ((positionX < X0) || (positionX > X0 +L0))
double posx (double t){
    double positionX = tmpposx;
    while(outOfBox){
        positionX -= signvc*L0;
    }
    return positionX;
}
void soild_fs(scalar fs, double tt){
    face vector fs_face[];
    vertex scalar phi[];
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
    fractions (phi, fs, fs_face);
//    fs.refine = fs.prolongation = fraction_refine;
    boundary({fs});
}

void mytracer (scalar f){
    vertex scalar phi[];
    face vector ff[];
    foreach_vertex() {
#ifdef MOVING_ON
        phi[] = (x < -3*rad && y < 0)? 1 : -1;
#else
//        phi[] = (x < -1.5*rad && y < 0)? 1 : -1;;
        phi[] = 0;
#endif
    }
    boundary ({phi});
    fractions (phi, f, ff);
}

int main(int argc, char * argv[])
{
    minlevel = 4;
    maxlevel = 9;
    eta_s=1e-2;
    eps = 1e-2;
    if (argc > 1) {
        maxlevel = atoi(argv[1]);
    }
    if (argc > 2) {
        eta_s = atof(argv[2]);
    }
    if (argc > 3) {
        eps = atof(argv[3]);
    }
    size (Ldomain);
#ifdef MOVING_ON
    origin (3.*rad - L0, -L0/2.);
#else
    origin (-3.*rad,     -L0/2.);
#endif
    DT = 1e-12;
    CFL=0.4;
    TOLERANCE = 1e-8;
    RELATIVE_RES_TOLERANCE = 0.1;
    NITERMAX=30;
    N = 1 << minlevel;
    const face vector muc[] = {2.0*rad/RE, 2.0*rad/RE, 2.0*rad/RE}; //Re=rho*U*2*rad/mu, mu=rho*U*2*rad/Re
    mu = muc;
    signvc = (vc.x > 0) ? 1 : (vc.x < 0)? -1 : 0;
    const vector U_sol[] = {vc.x, vc.y, vc.z};
    target_U = U_sol;

    fprintf(ferr, "INT_MAX= %d\n", INT_MAX);
    run();
}

scalar un[];

event init (t = 0) {
    if (!restore (file = "restart")) {
        int it = 0;
        do {
            soild_fs (fs, 0);
            mytracer (f);
        }while (adapt_wavelet({fs, f}, (double []){1e-5, 1e-5}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
        foreach() u.x[] = (1.0 - fs[]);
        boundary(all);
    }
    event("vtk_file");
}

event viscous_term (i++)
{
    soild_fs (fs, t+dt);
}

event set_dtmax (i++) {
    NITERMIN=1;
    NITERMAX=100;

    DT *= 1.05;
    DT = min(DT, CFL*Ldomain/pow(2, maxlevel+3));

    fprintf(ferr, "set_dtmax: tnext= %g Dt= %g", tnext, DT);
}


event end_timestep (i++){
    vorticity (u, omega);
    int curlevel=1;
    scalar l[]; foreach(reduction(max:curlevel)) {l[] = level; if (curlevel<level) curlevel=level;}
    #ifdef DEBUG_MODE
        vector gradu[], gradv[];
        gradients({u.x, u.y}, {gradu, gradv});
        foreach(){
            mod_du_dx[] = sqrt(sq(gradu.x[]) + sq(gradu.y[]) + sq(gradv.x[]) + sq(gradv.y[]));
        }
    #endif
    int tnc = count_cells();
    int nmax = (int)pow (2, maxlevel*dimension);
    fprintf (ferr, "i= %d t+dt= %g dt= %g eta_s= %g count_cells= %d nmax= %d compress_ratio= %g curlevel= %d\n", i, t+dt, dt, eta_s, tnc, nmax, (double)tnc/nmax, curlevel);
    double eps_arr[] = {1, 1, 1, 1};
    MinMaxValues((scalar *){p, u.x, u.y, omega}, eps_arr);
}

/**
We produce animations of the vorticity and tracer fields... */

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

//Output
//event vtk_file (i++){
event vtk_file (t += 0.01){
    char subname[80]; sprintf(subname, "v=0");
    scalar l[];
    vorticity (u, omega);
    foreach() {l[] = level; omega[] *= 1 - fs[]; }
//    vector mapped_data_lower[], mapped_data_upper[];
//    foreach() {
//        foreach_dimension()
//        {
//            mapped_data_lower.x[] = uf.x[];
//            mapped_data_upper.x[] = uf.x[1];
//        }
//    }
    #if DEBUG_BRINKMAN_PENALIZATION!=1
        output_vtu_MPI( (scalar *) {fs, f, omega, p, l, divutmpAfter}, (vector *) {u}, subname, 0 );
    #else
        output_vtu_MPI( (scalar *) {fs, f, omega, p, pf, l, divutmpAfter}, (vector *) {u, uf, dbp, total_rhs}, subname, 0 );
    #endif
}

//#define ADAPT_SCALARS {fs, f, u}
//#define ADAPT_EPS_SCALARS {1e-5, 1e-5, 1e-2, 1e-2}
event adapt (i++){
//    double eps_arr[] = ADAPT_EPS_SCALARS;
//    MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) {u.x, u.y, f, fs}, (double[]){eps, eps, eps, eps}, maxlevel = maxlevel, minlevel = minlevel);
}

event stop(t = 17);