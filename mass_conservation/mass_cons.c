#define BRINKMAN_PENALIZATION 1
#define DEBUG_BRINKMAN_PENALIZATION 1
#define DEBUG_OUTPUT_VTU_MPI
#define REDUCED 0
#define FILTERED
#include "../src_local/centered-weugene.h"
//#include "navier-stokes/double-projection.h"
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
int Nobst = 2; //250
scalar fs[], omega[];
double U0=0.01, rhol=1e+3, sig=73e-3, Lchar=5e-3, mul=1e-3, grav=-9.8;
double RE, CA, FR;//RE=500.0, CA=0.013, FR=20;
double Rrho=1000, Rmu=53.73, Ggrav;
double Radius_b = 0.125;
double obstacle_pattern(double x, double y, double xc, double yc, double size){
	return sq(x  - xc) + sq(y - yc) - sq(size);
}

void obstacles (scalar fs, int ns)
{
	face vector ffs[];
	double xc[ns], yc[ns], R[ns];
	srand (0);
	double dist = L0/ns;
	double size = 0.25*dist;
	for (int i = 0; i < ns; i++) {
		xc[i] = 0.25*L0;
		yc[i] = -0.5*L0 + 0.5*dist + dist*i;
		R[i] = size;
	}
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
						phi[] = intersection (phi[], obstacle_pattern(x, y, xc[i] - xp, yc[i] - yp, R[i]));
		phi[] = -phi[];
	}
	boundary ({phi});
	fractions (phi, fs, ffs);
}

void bubbles (scalar f)
{
	const int ns=1;
	face vector ff[];
	double xc[ns], yc[ns], R[ns];
	xc[0] = -0.250000; yc[0] = 0.000000; R[0] = Radius_b*L0;
//	xc[1] = -0.410100; yc[1] = -0.12890; R[1] = 0.02;
//	xc[2] = -0.378900; yc[2] =  0.14000; R[2] = 0.03;
//	xc[3] = -0.200000; yc[3] =  0.20000; R[3] = 0.10;
//	xc[4] =  0.000000; yc[4] =  0.26000; R[4] = 0.05;
//	xc[5] = -0.200000; yc[5] = -0.15000; R[5] = 0.05;

	vertex scalar phi[];
	foreach_vertex() {
		phi[] = HUGE;
		for (double xp = -L0; xp <= L0; xp += L0)
			for (double yp = -L0; yp <= L0; yp += L0)
				for (int i = 0; i < ns; i++)
					for (int i = 0; i < ns; i++)
						phi[] = intersection (phi[], (sq(x + xp - xc[i]) + sq(y + yp - yc[i]) - sq(R[i])));
		//phi[] = -phi[];
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
	eta_s = 1e-10;
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
			obstacles (fs, Nobst);
			bubbles(f);
			boundary (all);
		}while (adapt_wavelet({fs, f}, (double []){1e-4, 1e-4}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
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

//event snapshot (t += 0.5; t <= 10.8) {
//	char name[80];
//	sprintf (name, "snapshot-%g", t);
//	scalar pid[];
//	foreach()	pid[] = fmod(pid()*(npe() + 37), npe());
//	boundary ({pid});
//	dump (name);
//}
#undef SEPS
#define SEPS 1e-15
event logfile (i+=100)
{
	double avg = normf(u.x).avg, du = change (u.x, un)/(avg + SEPS);
	fprintf (ferr, "%d %d %g %g %d %d %d %d %d %d %.3g %.3g %.3g %.3g %.3g\n",
	maxlevel, i, t, dt,
	mgp.i, mgp.nrelax, mgp.minlevel,
	mgu.i, mgu.nrelax, mgu.minlevel,
	du, mgp.resa*dt, mgu.resa, statsf(u.x).sum, normf(p).max);
}

//Output
event vtk_file (t += 0.001; t<20){
	char subname[80]; sprintf(subname, "mc");
	scalar l[], npid[];
	vorticity (u, omega);
	foreach() {l[] = level; omega[] *= 1 - fs[]; npid[] = pid();}
	output_vtu_MPI( (scalar *) {fs, f, omega, p, l, npid. rho}, (vector *) {u, a}, subname, 1 );
}

#define ADAPT_SCALARS {f, fs, u.x, u.y}
#define ADAPT_EPS_SCALARS {1e-3, 1e-3, 1e-2, 1e-2}
event adapt (i++){
	double eps_arr[] = ADAPT_EPS_SCALARS;
	MinMaxValues(ADAPT_SCALARS, eps_arr);
	adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
	fs.refine = fs.prolongation = fraction_refine;
	boundary({fs});
}