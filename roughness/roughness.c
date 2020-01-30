/**
# Stokes flow through a complex porous medium

The medium is periodic and described using embedded boundaries.

This tests mainly the robustness of the representation of embedded
boundaries and the convergence of the viscous and Poisson
solvers. */
#define BRINKMAN_PENALIZATION 1
#define DEBUG_BRINKMAN_PENALIZATION 1
//#include "embed.h"
#include "../src_local/centered-weugene.h"
#include "two-phase.h"
#include "tension.h"
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

		for (double xp = -L0; xp <= L0; xp += L0)
		    phi[] = intersection (phi[], ( A*sin(2.0*pi*n*x/L0) + L0/pow(2, minlevel-1) - y ));
//		phi[] = -phi[];
	}
	boundary ({phi});

	fractions (phi, fs, ffs);
//	fractions_cleanup (fs, ffs);
}

/**
The domain is the periodic unit square centered on the origin. */

int main()
{
	origin (-0.5);
	eta_s = 1e-15;
	/**
	We turn off the advection term. The choice of the maximum timestep
	and of the tolerance on the Poisson and viscous solves is not
	trivial. This was adjusted by trial and error to minimize (possibly)
	splitting errors and optimize convergence speed. */

//	stokes = true;
	DT = 1e-3;

	N = 1 << maxlevel;
	rho1 = 1.; rho2 = 1./25.;
	mu1 = 1.; mu2 = 1.0/25.;
	f.sigma = 1.0e-2;
	for (scalar s in {f0,fs})
		s.refine = s.prolongation = fraction_refine;
	run();
}

scalar un[];

event init (t = 0) {
	if (!restore (file = "restart")) {
		int it = 0;
		do {
            rough_surface (fs, L0/20., 10);
			boundary (all); // this is necessary since BCs depend on embedded fractions
			foreach() {
				f[] = ( sq(x) + sq(y-L0/2.) < sq(0.3875*L0) ) ? 1 : 0;
				f[] = clamp(f[] + fs[], 0, 1);
			}
			boundary ({f});
		}while (adapt_wavelet({fs,f}, (double []){1e-3, 1e-3}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
	}
	/**
	We initialize the reference velocity. */
	foreach() un[] = u.x[];


}

/**
The gravity vector is aligned with the channel and viscosity is
unity. */
event acceleration (i++) {
	face vector av = a;
	foreach_face(y)	av.y[] -= 1;
}
/**
We check for a stationary solution. */

//event logfile (i++; i <= 500)
//{
//	double avg = normf_weugene(u.x,fs).avg, du = change_weugene (u.x, un, fs)/(avg + SEPS);
//	fprintf (ferr, "%d %d %d %d %d %d %d %d %.3g %.3g %.3g %.3g %.3g\n",
//	maxlevel, i,
//	mgp.i, mgp.nrelax, mgp.minlevel,
//	mgu.i, mgu.nrelax, mgu.minlevel,
//	du, mgp.resa*dt, mgu.resa, statsf_weugene(u.x).sum, normf(p).max);
//
//	/**
//	If the relative change of the velocity is small enough we stop this
//	simulation. */
//	if (i > 1 && (avg < 1e-9 || du < 1e-2)) {
//		/**
//		We are interested in the permeability $k$ of the medium, which is
//		defined by
//		$$
//		U = \frac{k}{\mu}\nabla p = \frac{k}{\mu}\rho g
//		$$
//		with $U$ the average fluid velocity.
//		*/
//		stats s = statsf_weugene (u.x, fs);
//		fprintf (fout, "%d %g\n", maxlevel, s.sum/s.volume);
//		fflush(fout);
//	}
//}

//event snapshot (t = 0.1; t += 0.01; t <= 10.8) {
//	char name[80];
//	sprintf (name, "snapshot-%g", t);
//	scalar pid[];
//	foreach()	pid[] = fmod(pid()*(npe() + 37), npe());
//					boundary ({pid});
//	dump (name);
//}

//Output
event vtk_file (i+=100;t<10){
	char subname[80]; sprintf(subname, "porous_BP");
	scalar l[];
	vorticity (u, omega);
	foreach() {l[] = level; omega[] *= 1 - fs[];}
	output_vtu_MPI( (scalar *) {fs, f, omega, p, l}, (vector *) {u,a}, subname, L0/pow(2., minlevel));
}

event adapt (i++) {
	adapt_wavelet ({f, fs, u}, (double[]){0.01, 0.01, 2e-3, 2e-3, 2e-3}, maxlevel=maxlevel, minlevel=minlevel);
}
/**
![Norm of the velocity field.](porous/nu-10.png)

![Pressure field.](porous/p-10.png)

![Adapted mesh, 10 levels of refinement.](porous/level-10.png)

~~~gnuplot Permeability as a function of resolution
set xlabel 'Level'
set grid
set ytics format '%.1e'
set logscale y
plot 'out' w lp t ''
~~~

~~~gnuplot Convergence history
set xlabel 'Iterations'
set logscale y
set ytics format '%.0e'
set yrange [1e-10:]
plot '../porous.ref' u 2:9 w l t '', '' u 2:10 w l t '', \
    '' u 2:11 w l t '', '' u 2:12 w l t '', '' u 2:13 w l t '', \
    'log' u 2:9 w p t 'du', '' u 2:10 w p t 'resp', \
    '' u 2:11 w p t 'resu', '' u 2:12 w p t 'u.x.sum', '' u 2:13 w p t 'p.max'
~~~

## See also

* [Stokes flow past a periodic array of cylinders](cylinders.c)
* [Stokes flow through a complex 3D porous medium](/src/examples/porous3D.c)
*/