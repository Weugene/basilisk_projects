/**
# Stokes flow through a complex porous medium

The medium is periodic and described using embedded boundaries.

This tests mainly the robustness of the representation of embedded
boundaries and the convergence of the viscous and Poisson
solvers. */
#define BRINKMAN_PENALIZATION 1
#define DEBUG_BRINKMAN_PENALIZATION 1
#define REDUCED 1
#define FILTERED
#include "../src_local/centered-weugene.h"
#include "two-phase.h"
#include "tension.h"
#if REDUCED
	#include "reduced.h"
#endif
#include "view.h"
#include "../src_local/output_vtu_foreach.h"
/**
We will vary the maximum level of refinement, starting from 5. */

int maxlevel = 9;
int minlevel = 5;
scalar f0[], fs[], omega[];
/**
The porous medium is defined by the union of a random collection of
disks. The number of disks can be varied to vary the porosity. */

void rough_surface (scalar fs, double A, int n)
{
	face vector ffs[];

	/**
	Once we have defined the random centers and radii, we can compute
	the levelset function $\phi$ representing the embedded boundary. */

	vertex scalar phi[];
	foreach_vertex() {
		phi[] = HUGE;

		/**
		Since the medium is periodic, we need to take into account all
		the disk images using periodic symmetries. */

//		for (double xp = -L0; xp <= L0; xp += L0)
		    phi[] = intersection (phi[], ( A*sin(2.0*pi*n*x/L0) + A + L0/pow(2, minlevel-2) - y ));
//		phi[] = -phi[];
	}
	boundary ({phi});

	fractions (phi, fs, ffs);
//	fractions_cleanup (fs, ffs);
}

void drop (scalar f, scalar fs)
{
	face vector ff[];
	vertex scalar phi[];
	foreach_vertex() {
		phi[] = HUGE;
		phi[] = intersection (phi[], (sq(x) + sq(y-L0/2.) - sq(0.25*L0) ));
		phi[] = -phi[];
	}
	boundary ({phi});
	fractions (phi, f, ff);
	foreach() f[] = clamp(f[] + fs[], 0, 1);
}

f[left] =neumann(0);
f[right]=neumann(0);
f[bottom] =dirichlet(1);
f[top]=neumann(0);

u.n[left] =neumann(0);
u.n[right]=neumann(0);
u.n[bottom] =dirichlet(0);
u.n[top]=neumann(0);

uf.t[left] =neumann(0);
uf.t[right]=neumann(0);
uf.t[bottom] =dirichlet(0);
uf.t[top]=neumann(0);

p[left] =neumann(0);
p[right]=neumann(0);
p[bottom] =neumann(0);
p[top]=dirichlet(0);

pf[left] =neumann(0);
pf[right]=neumann(0);
pf[bottom] =neumann(0);
pf[top]=dirichlet(0);

/**
The domain is the periodic unit square centered on the origin. */

int main()
{
    size(0.01);
	origin (-0.5*L0);
	eta_s = 1e-15;
	/**
	We turn off the advection term. The choice of the maximum timestep
	and of the tolerance on the Poisson and viscous solves is not
	trivial. This was adjusted by trial and error to minimize (possibly)
	splitting errors and optimize convergence speed. */

//	stokes = true;
	DT = 1e-4;
	TOLERANCE = 1e-6;
	//NITERMAX = 100;
	//mgp.nrelax = 100;
	N = 1 << 7;
	rho1 = 1000.; rho2 = 1.2;
	mu1 = 1.004e-3; mu2 = 18.5e-5;
	f.sigma = 73.e-3;
    #if REDUCED
	    G.y = -9.8;
	    Z.y = 0;
    #endif
	for (scalar s in {f0,fs})
		s.refine = s.prolongation = fraction_refine;
	run();
}

event init (t = 0) {
	if (!restore (file = "restart")) {
		int it = 0;
		do {
            rough_surface (fs, L0/80., 30);
			drop(f, fs);
			boundary (all); // this is necessary since BCs depend on embedded fractions
		}while (adapt_wavelet({fs,f}, (double []){1e-3, 1e-3}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
	}
	event("vtk_file");
}

/**
The gravity vector is aligned with the channel and viscosity is
unity. */
#if !REDUCED
event acceleration (i++) {
	face vector av = a;
	foreach_face(y)	av.y[] = -9.8;
}
#endif

event snapshot (t += 0.01; t <= 100) {
	char name[80];
	sprintf (name, "snapshot-%g", t);
	scalar pid[];
	foreach()	pid[] = fmod(pid()*(npe() + 37), npe());
					boundary ({pid});
	dump (name);
}

event logfile (i+=1)
{
    fprintf (ferr, "%d %d %g %g %d %d %d %d %d %d %.3g %.3g \n",
        maxlevel, i, t, dt,
        mgp.i, mgp.nrelax, mgp.minlevel,
        mgu.i, mgu.nrelax, mgu.minlevel,
        mgp.resa*dt, mgu.resa);
}

//Output
event vtk_file (t+=0.0001){
	char subname[80]; sprintf(subname, "rg");
	scalar l[];
	vorticity (u, omega);
	foreach() {l[] = level; omega[] *= 1 - fs[];}
	output_vtu_MPI( (scalar *) {fs, f, omega, p, l}, (vector *) {u, a}, subname, 0);
}

event adapt (i++) {
	adapt_wavelet ({f, fs, u}, (double[]){1e-3, 1e-3, 1e-2, 1e-2, 1e-2}, maxlevel=maxlevel, minlevel=minlevel);
    fs.refine = fs.prolongation = fraction_refine;
    boundary({fs});
}
