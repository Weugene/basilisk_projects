/**
# Stokes flow through a complex porous medium

The medium is periodic and described using embedded boundaries.

This tests mainly the robustness of the representation of embedded
boundaries and the convergence of the viscous and Poisson
solvers. */
//#define BRINKMAN_PENALIZATION 1
//#define DEBUG_BRINKMAN_PENALIZATION 1
#include "embed.h"
#include "../src_local/centered-weugene.h"
//#include "two-phase.h"
//#include "tension.h"
#include "view.h"
#include "../src_local/output_vtu_foreach.h"
/**
We will vary the maximum level of refinement, starting from 5. */

int maxlevel = 9;
int minlevel = 5;
scalar omega[];
/**
The porous medium is defined by the union of a random collection of
disks. The number of disks can be varied to vary the porosity. */

void porous (scalar cs, face vector fs)
{
	int ns = 250;
	double xc[ns], yc[ns], R[ns];
	srand (0);
	for (int i = 0; i < ns; i++)
		xc[i] = 0.5*noise(), yc[i] = 0.5*noise(), R[i] = 0.02 + 0.03*fabs(noise());

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
			for (double yp = -L0; yp <= L0; yp += L0)
				for (int i = 0; i < ns; i++)
				for (int i = 0; i < ns; i++)
					phi[] = intersection (phi[], (sq(x + xp - xc[i]) +
					                              sq(y + yp - yc[i]) - sq(R[i])));
		phi[] = -phi[];
	}
	boundary ({phi});

	fractions (phi, cs, fs);
	fractions_cleanup (cs, fs);
}

/**
The domain is the periodic unit square centered on the origin. */
/**
The boundary condition is zero velocity on the embedded boundary. */
u.n[embed] = dirichlet(0);
u.t[embed] = dirichlet(0);

int main()
{
	origin (-0.5, -0.5);
	periodic (right);
	periodic (top);

	/**
	We turn off the advection term. The choice of the maximum timestep
	and of the tolerance on the Poisson and viscous solves is not
	trivial. This was adjusted by trial and error to minimize (possibly)
	splitting errors and optimize convergence speed. */

	stokes = true;
	DT = 2e-5;
#if 1
	TOLERANCE = HUGE;
	NITERMIN = 2;
#else
	TOLERANCE = 1e-3;
    NITERMIN = 2;
#endif
	N = 1 << maxlevel;
	run();
}

scalar un[];

event init (t = 0) {
	if (!restore (file = "restart")) {
		int it = 0;
		do {
			porous (cs, fs);
			boundary (all); // this is necessary since BCs depend on embedded fractions
//			foreach() {
//				fs[] = sq(x - L0/2) + sq(y) < sq(Rsolid) ? 1 : 0;
//				f[] = f0[]*(x < length) + fs[];
//				u.x[] = f[];
//			}
//			boundary ({f,u.x});
		}while (adapt_wavelet({cs}, (double []){0.001}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
	}
	/**
	We initialize the reference velocity. */
	foreach() un[] = u.x[];


	/**
	The gravity vector is aligned with the channel and viscosity is
	unity. */
	const face vector g[] = {1.,0.};
	a = g;
	mu = fm;
}

/**
We check for a stationary solution. */

event logfile (i++; i <= 500)
{
	double avg = normf(u.x).avg, du = change (u.x, un)/(avg + SEPS);
	fprintf (ferr, "%d %d %d %d %d %d %d %d %.3g %.3g %.3g %.3g %.3g\n",
	maxlevel, i,
	mgp.i, mgp.nrelax, mgp.minlevel,
	mgu.i, mgu.nrelax, mgu.minlevel,
	du, mgp.resa*dt, mgu.resa, statsf(u.x).sum, normf(p).max);

	/**
	If the relative change of the velocity is small enough we stop this
	simulation. */
	if (i > 1 && (avg < 1e-9 || du < 1e-2)) {
		/**
		We are interested in the permeability $k$ of the medium, which is
		defined by
		$$
		U = \frac{k}{\mu}\nabla p = \frac{k}{\mu}\rho g
		$$
		with $U$ the average fluid velocity.
		*/

		stats s = statsf (u.x);
		fprintf (fout, "%d %g\n", maxlevel, s.sum/s.volume);
		fflush(fout);

		/**
		We output fields and dump the simulation. */
		scalar nu[];
		foreach() nu[] = sqrt (sq(u.x[]) + sq(u.y[]));
		boundary ({nu});

		view (fov = 19.3677);

		draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
		squares ("nu", linear = true, spread = 8);
		char name[80];
		sprintf (name, "nu-%d.png", maxlevel);
		save (name);

		draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
		squares ("p", linear = false, spread = -1);
		sprintf (name, "p-%d.png", maxlevel);
		save (name);

		draw_vof ("cs", "fs", filled = -1, fc = {1,1,1});
		squares ("level");
		sprintf (name, "level-%d.png", maxlevel);
		save (name);

		sprintf (name, "dump-%d", maxlevel);
//		dump (name);

		/**
		We stop at level 10. */

		if (maxlevel == 10) return 1; /* stop */

		/**
		We refine the converged solution to get the initial guess for the
		finer level. We also reset the embedded fractions to avoid
		interpolation errors on the geometry. */

		maxlevel++;
		#if 0
			refine (level < maxlevel && cs[] > 0. && cs[] < 1.);
		#else
			adapt_wavelet ({cs,u}, (double[]){1e-2,2e-6,2e-6}, maxlevel);
		#endif
		porous (cs, fs);
		boundary (all); // this is necessary since BCs depend on embedded fractions
	}
}

//Output
event vtk_file (i+=1){
	char subname[80]; sprintf(subname, "porous_BP");
	scalar l[];
	vorticity (u, omega);
	foreach() {l[] = level; omega[] *= 1 - cs[];}
	output_vtu_MPI( (scalar *) {cs, omega, p, l}, (vector *) {u}, subname, L0/pow(2., minlevel));
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