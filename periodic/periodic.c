#define REDUCED 1
#define FILTERED
#include "navier-stokes/centered.h"
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
#include "tension.h"
#if REDUCED
	#include "reduced.h"
#endif

#include "view.h"
#include "../src_local/output_vtu_foreach.h"

int maxlevel = 8;
int minlevel = 5;
scalar omega[];
double U0=0.01, rhol=1e+3, sig=73e-3, Lchar=5e-3, mul=1e-3, grav=-9.8;
double RE, CA, FR;//RE=500.0, CA=0.013, FR=20;
double Rrho=1000, Rmu=53.73, Ggrav;
double Radius_b = 0.125;

void bubbles (scalar f)
{
	const int ns=1;
	face vector ff[];
	double xc[ns], yc[ns], R[ns];
	xc[0] = -0.250000; yc[0] = 0.000000; R[0] = Radius_b*L0;

	vertex scalar phi[];
	foreach_vertex() {
		phi[] = HUGE;
		for (double xp = -L0; xp <= L0; xp += L0)
			for (double yp = -L0; yp <= L0; yp += L0)
				for (int i = 0; i < ns; i++)
					for (int i = 0; i < ns; i++)
						phi[] = intersection (phi[], (sq(x + xp - xc[i]) + sq(y + yp - yc[i]) - sq(R[i])));
	}
	boundary ({phi});
	fractions (phi, f, ff);
}
/**
The domain is the periodic unit square centered on the origin. */


int main(int argc, char * argv[])
{
    if (argc > 1) {
        Radius_b = atof(argv[1]); //convert from string to float
    }
    if (argc > 2) {
        maxlevel = atoi(argv[2]); //convert from string to float
    }
    size (1.0);
	origin (-0.5*L0, -0.5*L0);
	periodic (right);
	periodic (top);
	DT = 1e-3;
	TOLERANCE = 1e-8;
	N = 1 << maxlevel;
	RE=U0*Lchar*rhol/mul; CA=U0*mul/sig; FR=sq(U0)/(grav*Lchar);
	rho1 = 1.; rho2 = rho1/Rrho;
	mu1 = 1./RE; mu2 = mu1/Rmu;
	f.sigma = 1./RE/CA;
	Ggrav = 1./FR;
#if REDUCED
	G.x = Ggrav;
	Z.x = 0;
#endif
	fprintf(ferr, "RE=%g CA=%g FR=%g \n"
			   "mu1=%g mu2=%g rho1=%g rho2=%g \n"
	           "sigma=%g grav=%g L0=%g\n",
	           RE, CA, FR, mu1, mu2, rho1, rho2, f.sigma, Ggrav, L0);
	run();
}

scalar un[];

event init (t = 0) {
	if (!restore (file = "restart")) {
		int it = 0;
		do {
			bubbles(f);
			boundary (all);
		}while (adapt_wavelet({f}, (double []){1e-4, 1e-4}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
	}
	event("vtk_file");
}

/**
The gravity vector is aligned with the channel and viscosity is
unity. */

#if !REDUCED
event acceleration (i++) {
	face vector av = a;
	foreach_face(x)	av.x[] += Ggrav;
}
#endif

event logfile (i+=100)
{
    double avgf = normf(f).avg;
	fprintf (ferr, "%d %d %g %g %g\n",
	maxlevel, i, t, dt, sq(L0) - avgf);
}


//Output
//event vtk_file (i++; t<20){
event vtk_file (t += 0.01; t<20){
	char subname[80]; sprintf(subname, "mc");
	scalar l[], npid[];
	vorticity (u, omega);
	foreach() {l[] = level;}
	output_vtu_MPI( (scalar *) {f, omega, p, l}, (vector *) {u, a}, subname, 1 );
}

#define ADAPT_SCALARS {f, u.x, u.y}
#define ADAPT_EPS_SCALARS {1e-3, 1e-2, 1e-2}
event adapt (i++){
	double eps_arr[] = ADAPT_EPS_SCALARS;
//	MinMaxValues(ADAPT_SCALARS, eps_arr);
	adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);

}
//plot 'plot' using 3:(a=max(a,$5),0/0) notitle, 'plot' using 3:($5/a) with linespoints