#define BRINKMAN_PENALIZATION 1
#define DEBUG_BRINKMAN_PENALIZATION 1
#define DEBUG_OUTPUT_VTU_MPI
#define REDUCED 0
#define FILTERED
#include "../src_local/centered-weugene.h"
#include "view.h"
#include "../src_local/output_vtu_foreach.h"
#include "tracer.h"

scalar f[];
scalar * tracers = {f};
face vector muv[];

int maxlevel = 10;
int minlevel = 4;
scalar fs[], omega[];

/**
The domain is the periodic unit square centered on the origin. */

/**
The fluid is injected on the left boundary with a unit velocity. The
tracer is injected in the lower-half of the left boundary. An outflow
condition is used on the right boundary. */

u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(y < 0);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

int main(int argc, char * argv[])
{
    if (argc > 1) {
        maxlevel = atoi(argv[2]); //convert from string to float
    }
    size (8.0);
    origin (-1.5, -L0/2.);
	eta_s = 1e-10;
	TOLERANCE = 1e-8;
    N = 512;
    mu = muv;
    run();
}

scalar un[];

event init (t = 0) {
    /**
    The domain is the intersection of a channel of width unity and a
    circle of diameter 0.125. */
    fraction (fs, sq(0.5) - sq(x) - sq(y));

    /**
    We set the initial velocity field. */

    foreach() u.x[] = 1 - fs[];

//
//	if (!restore (file = "restart")) {
//		int it = 0;
//		do {
//			fraction (fs, sq(x - 0.5) + sq(y) - sq(0.25));
//		}while (adapt_wavelet({fs, f}, (double []){1e-4, 1e-4}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
//	}
//	event("vtk_file");
}

/**
We set a constant viscosity corresponding to a Reynolds number of 160,
based on the cylinder diameter (0.125) and the inflow velocity (1). */

event properties (i++) {
    foreach_face() muv.x[] = fm.x[]*0.125/160.;
}

event logfile (i++) {
    fprintf (ferr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
}

/**
We produce animations of the vorticity and tracer fields... */

event movies (i += 4; t <= 15.) {
    scalar omega[], m[];
    vorticity (u, omega);
    foreach() m[] = 0.5 - fs[]; boundary ({m});
    output_ppm (omega, file = "vort.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
            min = -10, max = 10, linear = true, mask = m);
    output_ppm (f, file = "f.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
            linear = false, min = 0, max = 1, mask = m);
}

//Output
event vtk_file (t += 0.1){
	char subname[80]; sprintf(subname, "rk");
	scalar l[];
	vorticity (u, omega);
	foreach() {l[] = level; omega[] *= 1 - fs[]; }

#if BRINKMAN_PENALIZATION==1
	output_vtu_MPI( (scalar *) {fs, f, omega, p, l}, (vector *) {u}, subname, 0 );
#else
    output_vtu_MPI( (scalar *) {fs, f, omega, p, l}, (vector *) {u, target_U, dbp, total_rhs, utau, grad_utau_n}, subname, 0 );
#endif
}

#define ADAPT_SCALARS {fs, f, u}
#define ADAPT_EPS_SCALARS {1e-2, 3e-2, 3e-2, 3e-2}
event adapt (i++){
	double eps_arr[] = ADAPT_EPS_SCALARS;
//	MinMaxValues(ADAPT_SCALARS, eps_arr);
	adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
	fs.refine = fs.prolongation = fraction_refine;
	boundary({fs});
}

