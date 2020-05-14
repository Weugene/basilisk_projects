#define BRINKMAN_PENALIZATION 1
#define DEBUG_BRINKMAN_PENALIZATION 1
#define DEBUG_OUTPUT_VTU_MPI
#define FILTERED
#define RELATIVE_RESIDUAL
#define EPS_MAXA 2
#define MODIFIED_CHORIN 1
face vector fs_face[];
(const) face vector target_Uf = zerof;
#include "../src_local/centered-weugene.h"
//#include "../src_local/three-phase-weugene.h"//rho3, mu3
#include "view.h"
#include "../src_local/output_vtu_foreach.h"
#include "tracer.h"
scalar f[];
scalar * tracers = {f};

int maxlevel = 11;
int minlevel = 4;
double xx0 = 0, rad = 0.00625, RE=40.;
scalar fs[], omega[], divu[];
double signvc = -1;
coord vc = {-1.0, 0.0, 0.0};

/**
The domain is the periodic unit square centered on the origin. */

/**
The fluid is injected on the left boundary with a unit velocity. The
tracer is injected in the lower-half of the left boundary. An outflow
condition is used on the right boundary. */

u.n[left]  = neumann(0.);
p[left]    = dirichlet(0.);
pf[left]   = dirichlet(0.);
f[left]   = dirichlet(y < 0);

u.n[right] = neumann(0);
p[right]   = neumann(0.);
pf[right]  = neumann(0.);

#define tmpposx (xx0 + vc.x*t)
#define outOfBox ((positionX < X0) || (positionX > X0 +L0))
double posx (double t){
    double positionX = tmpposx;
    while(outOfBox){
        positionX -= signvc*L0;
    }
    return positionX;
}
void soild_fs(scalar fs, face vector fs_face, double tt){
    double x00 = posx(tt);

    vertex scalar phi[];
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
    fs.refine = fs.prolongation = fraction_refine;
    boundary({fs});
}

void mytracer (scalar f){
    vertex scalar phi[];
    face vector ff[];
    foreach_vertex() {
        phi[] = (x < -3*rad && y < 0)? 1 : -1;
    }
    boundary ({phi});
    fractions (phi, f, ff);
}

int main(int argc, char * argv[])
{
    if (argc > 1) {
        maxlevel = atoi(argv[1]); //convert from string to float
    }
    size (1.0);
    origin (3.*rad - L0, -L0/2.);
    eta_s = 1e-6;
    DT = 1e-8;
    CFL=0.4;
    TOLERANCE = 1e-8;
    RELATIVE_RES_TOLERANCE = 0.1;
    NITERMAX=30;
    N = 512;
    const face vector muc[] = {2.0*rad/RE, 2.0*rad/RE, 2.0*rad/RE}; //Re=rho*U*2*rad/mu, mu=rho*U*2*rad/Re
    mu = muc;
    signvc = (vc.x > 0) ? 1 : (vc.x < 0)? -1 : 0;
    const vector U_sol[] = {vc.x, vc.y, vc.z};
    target_U = U_sol;
    const face vector U_solf[] = {vc.x, vc.y, vc.z};
    target_Uf = U_solf;
    fprintf(ferr, "INT_MAX= %d\n", INT_MAX);
    run();
}

scalar un[];

event init (t = 0) {
    if (!restore (file = "restart")) {
        int it = 0;
        do {
            soild_fs (fs, fs_face, 0);
            mytracer (f);
        }while (adapt_wavelet({fs, f}, (double []){1e-5, 1e-5}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
        foreach() u.x[]= 0;
        boundary(all);
    }
}

event set_dtmax (i++) {
    if (i<=100) {
        NITERMIN=100;
        NITERMAX=150;
    }else{
        NITERMIN=10;
        NITERMAX=30;
    }
    DT *= 1.05;
    DT = min(DT, CFL/pow(2, maxlevel+3));
    fprintf(ferr, "set_dtmax: tnext= %g Dt= %g", tnext, DT);
}

event vof(i++){
    soild_fs (fs, fs_face, t + 0.5*dt);
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
    double avg = normf(u.x).avg;
    fprintf (ferr, "i=%d t=%g dt=%g i_p=%d i_u=%d divu_max=%g avg= %g\n", i, t, dt, mgp.i, mgu.i, Linfu, avg);
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
    char subname[80]; sprintf(subname, "rk");
    scalar l[];
    vorticity (u, omega);
    foreach() {l[] = level; omega[] *= 1 - fs[]; }
    vector mapped_data_lower[], mapped_data_upper[];
    foreach() {
        foreach_dimension()
        {
            mapped_data_lower.x[] = uf.x[];
            mapped_data_upper.x[] = uf.x[1];
        }
    }
    #if DEBUG_BRINKMAN_PENALIZATION!=1
    output_vtu_MPI( (scalar *) {fs, f, omega, p, l, divu}, (vector *) {u, mapped_data_lower, mapped_data_upper}, subname, 0 );
    #else
    output_vtu_MPI( (scalar *) {fs, f, omega, p, l, divu}, (vector *) {u, dbp, total_rhs, mapped_data_lower, mapped_data_upper}, subname, 0 );
    #endif
}

#define ADAPT_SCALARS {fs, f, u}
#define ADAPT_EPS_SCALARS {1e-5, 1e-5, 1e-2, 1e-2}
event adapt (i++){
    double eps_arr[] = ADAPT_EPS_SCALARS;
    MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
}
event stop(t = 17);