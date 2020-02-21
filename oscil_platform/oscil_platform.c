#define DEBUG_MINMAXVALUES
#define DEBUG_OUTPUT_VTU_MPI
#define REDUCED 1
#define FILTERED
#include "../src_local/centered-weugene.h"
#include "contact.h"
#include "vof.h"
#include "two-phase.h"
#include "tension.h"
#if REDUCED
	#include "reduced.h"
#endif
#include "../src_local/output_vtu_foreach.h"


int maxlevel = 9;
int minlevel = 4;
double A=15.0, freq=100, Vmax;
double diam=2e-2;
double rhoL1=1000, rhoL2=1.2;
double muL1=1e-3, muL2=3.124e-5;
double SIGMA=73e-3, grav=9.8066;
double lcap, RE, BO;

static double k_a_wave[6][2] = {
		{28e+3,  4.375},
		{32.5e+3,3.777},
		{35e+3,  3.960},
		{48e+3,  12.506},
		{60.9e+3,19.760},
		{85e+3,  41.953}
};
/**
To set the contact angle, we allocate a [height-function
field](/src/heights.h) and set the contact angle boundary condition on
its tangential component. */

vector h[];
double theta0 = 30;
h.t[bottom] = contact_angle (theta0*pi/180.);

u.n[top]  = neumann(0);
u.t[top]  = neumann(0);
uf.n[top]  = neumann(0);
uf.t[top]  = neumann(0);
p[top]    = dirichlet(0);
pf[top]    = dirichlet(0);
f[top]    = (u.n[] >= 0) ? neumann(0) : 0;

u.n[bottom] = dirichlet(0);
u.t[bottom] = dirichlet(0);
uf.n[bottom] = dirichlet(0);
uf.t[bottom] = dirichlet(0);
p[bottom]   = neumann(0);
pf[bottom]   = neumann(0);
f[bottom]    = neumann(0);

u.n[left] = dirichlet(0);
u.t[left] = dirichlet(0);//dirichlet(Vmax*(1. - cos(2*pi*freq*t)));
uf.n[left] = dirichlet(0);
uf.t[left] = dirichlet(0);
p[left]   = neumann(0);
pf[left]   = neumann(0);
f[left]    = neumann(0);

u.n[right] = dirichlet(0);
u.t[right] = dirichlet(0);//dirichlet(Vmax*(1. - cos(2*pi*freq*t)));
uf.n[right] = dirichlet(0);
uf.t[right] = dirichlet(0);
p[right]   = neumann(0);
pf[right]   = neumann(0);
f[right]    = neumann(0);

int main()
{
	/**
	We use a constant viscosity and density. */
	mu1 = muL1; mu2 = muL2; rho1 = rhoL1; rho2 = rhoL2;
	/**
	We must associate the height function field with the VOF tracer, so
	that it is used by the relevant functions (curvature calculation in
	particular). */
	f.height = h;
	/**
	We set the surface tension coefficient and run for the range of
	contact angles. */
	DT=1e-5;
	f.sigma = SIGMA;
	Vmax = A*grav/(2*pi*freq);
	lcap = sqrt(SIGMA/(fabs(rho1-rho2)*grav));
	size (5.*lcap);//0.231e-3 m
	fprintf(ferr, "lc=%g Re1min=%g Re1max=%g Re2min=%g Re2max=%g Bomin=%g Bomax=%g", lcap,
	        rho1*freq/(mu1*sq(k_a_wave[5][0])), rho1*freq/(mu1*sq(k_a_wave[0][0])),
	        rho2*freq/(mu2*sq(k_a_wave[5][0])), rho2*freq/(mu2*sq(k_a_wave[0][0])),
	        1./sq(k_a_wave[5][0]*lcap),         1./sq(k_a_wave[0][0]*lcap)
			);
#if REDUCED
	G.y = 0;
	Z.y = 0;
#endif
	run();
}

/**
The initial drop is a quarter of a circle. */

event init (t = 0)
{
	//scalar fbubble[];
	if (!restore (file = "restart")) {
		int it = 0;
		do {
			fraction (f, -y + 0.25*L0 + 0.01*L0*sin(2.*pi*x/(0.1*L0)));
			//!fraction (fbubble, sq(0.5*diam) - sq(x-L0/2.) -sq(y-0.005));
			boundary (all); // this is necessary since BCs depend on embedded fractions
		}while (adapt_wavelet({f}, (double []){1e-4}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
	}
}

event acceleration (i++) {
#if REDUCED
	G.y = A*grav*sin(2.0*pi*freq*t);//m^2/s
	Z.y = 0;
#endif
}

#if 1
event logfile (i++)
{
  fprintf (fout, "%g %g\n", t, normf(u.x).max);
}

event snapshot (t +=10./freq) {
	char name[80];
	sprintf (name, "snapshot-%g", t);
	scalar pid[];
	foreach()	pid[] = fmod(pid()*(npe() + 37), npe());
			boundary ({pid});
	dump (name);
}
#endif

/**
At equilibrium (t = 10 seems sufficient), we output the interface
shape and compute the (constant) curvature. */

event end (t = 10)
{
	output_facets (f, stdout);

	scalar kappa[];
	curvature (f, kappa);
	stats s = statsf (kappa);
	double R = s.volume/s.sum, V = 2.*statsf(f).sum;
	fprintf (ferr, "%d %g %.5g %.3g\n", N, theta0, R/sqrt(V/pi), s.stddev);
}

scalar omega[];
event vtk_file (t+=1./(20.*freq)){
	char subname[80]; sprintf(subname, "osc");
	scalar l[];
	vector aa[];
	vorticity (u, omega);
	foreach() {
		l[] = level;
		foreach_dimension() aa.x[] = G.x;
	}
	output_vtu_MPI( (scalar *) {f, omega, p, l}, (vector *) {u, aa}, subname, 0);
}

#define ADAPT_SCALARS {f, u.x, u.y}
#define ADAPT_EPS_SCALARS {1e-3, 1e-2, 1e-2}
event adapt (i++){
	double eps_arr[] = ADAPT_EPS_SCALARS;
    MinMaxValues(ADAPT_SCALARS, eps_arr);
	adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
}