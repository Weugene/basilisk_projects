#define BRINKMAN_PENALIZATION 1
#define DEBUG_BRINKMAN_PENALIZATION 1
#define DEBUG_OUTPUT_VTU_MPI
//#define FILTERED
#include "../src_local/centered-weugene.h"
#include "view.h"
#include "../src_local/output_vtu_foreach.h"
#include "tracer.h"
scalar f[];
scalar * tracers = {f};

int maxlevel = 12;
int minlevel = 4;
double xx0 = 7, rad = 0.0625;
scalar fs[], omega[], divu[];

/**
The domain is the periodic unit square centered on the origin. */

/**
The fluid is injected on the left boundary with a unit velocity. The
tracer is injected in the lower-half of the left boundary. An outflow
condition is used on the right boundary. */

u.n[left]  = neumann(0.);
p[left]    = dirichlet(0.);
pf[left]   = dirichlet(0.);


u.n[right] = dirichlet(-1.);
p[right]   = neumann(0.);
pf[right]  = neumann(0.);
f[right]    = dirichlet(y < 0);
void soild_fs(scalar fs, double t){
    fraction (fs, sq(rad) - sq(x - xx0) - sq(y));
    fs.refine = fs.prolongation = fraction_refine;
    boundary({fs});
}

int main(int argc, char * argv[])
{
	if (argc > 1) {
		maxlevel = atoi(argv[2]); //convert from string to float
	}
	size (8.0);
	origin (-0.5, -L0/2.);
	eta_s = 1e-10;
    DT=1e-2;
	TOLERANCE = 1e-8;
	N = 512;
    const face vector muc[] = {0.00078125,0.00078125};
    mu = muc;
	run();
}

scalar un[];

event init (t = 0) {
	if (!restore (file = "restart")) {
		int it = 0;
		do {
            soild_fs (fs, 0);
		}while (adapt_wavelet({fs, f}, (double []){1e-3, 1e-3}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
	    boundary(all);
	}
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
    fprintf (ferr, "i=%d t=%g dt=%g iter_p=%d iter_u=%d div u=%g \n", i, t, dt, mgp.i, mgu.i, Linfu);
}

/**
We produce animations of the vorticity and tracer fields... */

event images (t += 0.1) {
    static FILE * fp = popen ("ppm2gif > vort.gif", "w");
    vorticity (u, omega);
    /**
    Cells for which *m* is negative will be black in the movie. */
    scalar m[];
    foreach()
            m[] = 0.5 - fs[];
                    boundary ({m});
    output_ppm (omega, fp, box = {{-0.5,-0.5},{7.5,0.5}}, mask = m,
            min=-10, max=10, linear=true);
}

//Output
event vtk_file (t += 0.01){
    char subname[80]; sprintf(subname, "rk");
    scalar l[];
    vorticity (u, omega);
    foreach() {l[] = level; omega[] *= 1 - fs[]; }

    #if DEBUG_BRINKMAN_PENALIZATION!=1
        output_vtu_MPI( (scalar *) {fs, f, omega, p, l, divu}, (vector *) {u}, subname, 0 );
    #else
        output_vtu_MPI( (scalar *) {fs, f, omega, p, l, divu}, (vector *) {u, dbp, total_rhs}, subname, 0 );
    #endif
}

#define ADAPT_SCALARS {fs, f, u}
#define ADAPT_EPS_SCALARS {1e-3, 1e-3, 3e-2, 3e-2}
event adapt (i++){
    double eps_arr[] = ADAPT_EPS_SCALARS;
    //	MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
    fs.refine = fs.prolongation = fraction_refine;
    boundary({fs});
}
event stop(t = 7);