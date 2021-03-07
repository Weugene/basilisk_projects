#define BRINKMAN_PENALIZATION 1
#define DEBUG_BRINKMAN_PENALIZATION 1
#define DEBUG_MINMAXVALUES
#define DEBUG_OUTPUT_VTU_MPI
#define REDUCED 1
#define FILTERED
#include "../src_local/centered-weugene.h"
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
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
int minlevel = 4;
#define diam 5e-3//meter 5e-3
double Q = 1.e-5;//cm^3/s Q=v*S
#define v_in 4.0*Q/(pi*sq(diam))//inlet velocity
double hjet = 3e-2;//meter
double length = 2.55e-2;//4.1e-3
double hbreak = 1.0e-2, hhump=5e-3;
double Tac = 1e-1;

double SIGMA = 30e-3;//30-70e-3
double MUL = 3.5e-3, MUG = 18.5e-4;
double RHOL = 1e+3, RHOG = 1.2;
double grav =-9.8;
scalar f0[], fs[], omega[];

u.n[top]  = f0[] > 0 ? dirichlet( - (sqrt(sq(v_in) + 0*2.*grav*length))*f0[]) : neumann(0);
//u.n[top]  = f0[] > 0 ? dirichlet( - v_in*f0[]*7(t/Tac, 1)) : neumann(0);
u.t[top]  = f0[] > 0 ? dirichlet(0) : neumann(0);
uf.n[top]  = f0[] > 0 ? dirichlet( - (sqrt(sq(v_in) + 0*2.*grav*length))*f0[]) : neumann(0);
//uf.n[top]  = f0[] > 0 ? dirichlet( - v_in*f0[]*min(t/Tac, 1)) : neumann(0);
uf.t[top]  = f0[] > 0 ? dirichlet(0) : neumann(0);
//u.n[top]  = dirichlet( - v_in*f0[]*min(t/Tac, 1));
//u.t[top]  = dirichlet(0);
//uf.n[top]  = dirichlet( - v_in*f0[]*min(t/Tac, 1));
//uf.t[top]  = dirichlet(0);
p[top]    = neumann(0);
pf[top]    = neumann(0);
f[top]    = (f0[] > 0) ? dirichlet(f0[]) : (u.n[] >= 0) ? neumann(0) : 0;

u.n[bottom] = dirichlet(0);
u.t[bottom] = dirichlet(0);
uf.n[bottom] = dirichlet(0);
uf.t[bottom] = dirichlet(0);
p[bottom]   = dirichlet(0);
pf[bottom]   = dirichlet(0);
f[bottom]    = neumann(0);

u.n[left] = neumann(0);
u.t[left] = neumann(0);
uf.n[left] = neumann(0);
uf.t[left] = neumann(0);
p[left]   = dirichlet(0);
pf[left]   = dirichlet(0);
f[left]    = ((y < - hbreak) || (u.n[] <= 0)) ? neumann(0) : 0;

u.n[right] = neumann(0);
u.t[right] = neumann(0);
uf.n[right] = neumann(0);
uf.t[right] = neumann(0);
p[right]   = dirichlet(0);
pf[right]   = dirichlet(0);
f[right]    =  ((y < - hbreak) || (u.n[] >= 0)) ? neumann(0) : 0;
/**
The porous medium is defined by the union of a random collection of
disks. The number of disks can be varied to vary the porosity. */
//double f_surf(double x){
//	return 0.5*(L0 - 2*hjet)*(-1 - tanh(100*(sq(x) - sq(0.4*L0))/sq(L0)));
////	return 1 - pow((y+L0)/L0,n) - pow(x/(0.4*L0),n);
////	return 1 - pow((y+L0)/L0,n) - pow(x/(0.4*L0),n);
////	return pow(pow(R,n) -pow(x,n), 1./n) - R;
//}
double f_surf(double x, double y, double R, int n){
	return 0.5*hbreak*(-1.0 - tanh(100.0*(sq(x) - sq(0.37*L0))/sq(L0)));
//	return pow(pow(R,n) -pow(x,n), 1./n) - R;
//	return ((fabs(x)<0.4*L0 && y<=0) ? 1 : (y<=5.*fabs(x)-2.0*L0) ? 1 : -1);
}

double f_surf_hump(double x, double y){
	return 0.5*hbreak*(-1.0 - tanh(100.0*(sq(x) - sq(0.37*L0))/sq(L0))) + hhump*exp(-10000.*sq(fabs(x)-0.35*L0)/sq(L0));
//	return pow(pow(R,n) -pow(x,n), 1./n) - R;
//	return ((fabs(x)<0.4*L0 && y<=0) ? 1 : (y<=5.*fabs(x)-2.0*L0) ? 1 : -1);
}
void substrate (scalar fs)
{
//	face vector ffs[];
	vertex scalar phi[];
	foreach_vertex()
//	phi[] = (f_surf(x, y, 0, 0) > y)? 1 : 0;
	phi[] = (f_surf_hump(x, y) > y)? 1 : 0;
	//phi[] = ((fabs(x)<0.4*L0 && y<=0) ? 1 : (fabs(x)>0.4*L0 && y<=-5.*fabs(x)+2.0*L0) ? 1 : 0);
//	boundary ({phi});
	foreach() fs[] = 0.25*(phi[] + phi[1,0] + phi[0,1] + phi[1,1]);
//	fraction (fs, (f_surf(x, y, 0.26, 10)));
//	fraction (fs, (  4*     -0.5*(L0 - 2.0*hjet)*(1.0 + tanh(50.0*(x*x - 0.37*L0*0.37*L0)/(L0*L0))) - y));
//	foreach() fprintf(ferr, "fs=%g", fs[]);
//	fraction (fs, y - f_surf(x));
}

double f_jet(double x, double y, int n){
//	if (y>hjet-length && fabs(x)<0.5*diam) fprintf(ferr, "%g <? %g  %g >? %g\n", fabs(x), 0.5*diam*sqrt(v_in/sqrt(sq(v_in) + 2.*fabs(grav)*fabs(hjet-y))), y, hjet-length);
	return 0.5*diam*sqrt(v_in/sqrt(sq(v_in) + 2.*fabs(grav)*fabs(hjet-y)));
	//return hjet - length*pow(1 -pow(x/(0.5*diam),n), 1./n);
}

void jet_in (scalar f,scalar fs)
{
	double Rsurf = 0.34*L0;
	double R = (sq(diam) + sq(Rsurf))/(2.0*diam);
	double xc = 0, yc = -R + diam;
	vertex scalar phi[];
	foreach_vertex()
	phi[] = ((sq(R) - sq(x - xc) - sq(y - yc) >= 0 && fabs(x) <= Rsurf) ||
	         (fabs(x) <= f_jet(x, y, 0) && y >= hjet - length) ||
	         (y <= - 0.5*hbreak ) ) ? 1 : 0;
	boundary ({phi});
	foreach() {
		f0[] = 0.25*(phi[] + phi[1,0] + phi[0,1] + phi[1,1]);
		f[]=clamp(f0[] + fs[], 0, 1);
	}

//	fraction (f0, sq(0.5*diam) - sq(x));//sq(0.5*diam) - sq(y - hjet) - sq(x)
//	foreach() f1[] = f0[]*(y > hjet-length);//0.5*hjet
//	fraction (f0, 1 - pow((y-hjet)/length,10) - pow(x/(0.5*diam),10) );//sq(0.5*diam) - sq(y - hjet) - sq(x)
//	fraction (fc, (sq(R) - sq(x - xc) - sq(y - yc) > 0 && y > 0)? 1 : -1 );
//	foreach() {f0[] =clamp(f0[] + fc[], 0, 1); f1[] = f0[];} //clamp(f0[]+fs1[], 0, 1);
}

int main(int argc, char * argv[])
{
	if (argc > 1) {
		Q = atof(argv[1]); //convert from string to float
	}
	if (argc > 2) {
		maxlevel = atoi(argv[2]); //convert from string to float
	}
	L0 = 0.16; //meter
	origin (-0.5*L0, -(L0-hjet) );
	eta_s = 1e-15;
	DT = 2e-4;
	TOLERANCE = 1e-6;
//	NITERMAX = 100;
	N = 1 << 8;
	rho1 = RHOL; rho2 = RHOG;
	mu1 = MUL; mu2 = MUG;
	f.sigma = SIGMA;
#if REDUCED
	G.y = grav;
	Z.y = 0;
#endif
	run();
}

event init (t = 0) {
	if (!restore (file = "restart")) {
		int it = 0;
		do {
			jet_in(f, fs);
			foreach() {
				u.y[] = - v_in*f[]*(y - hjet + length - 0.01*L0> 0.0 );
			}
			f.refine = f.prolongation = fraction_refine;
			fs.refine = fs.prolongation = fraction_refine;
			boundary (all); // this is necessary since BCs depend on embedded fractions
		}while (adapt_wavelet({fs,f}, (double []){1e-5, 1e-4}, maxlevel=maxlevel+1, minlevel=minlevel).nf != 0 && ++it <= 10);
		refine( fabs(x) < 0.8*diam && fabs(y-hjet) < 1.1*length && level < maxlevel );
		refine( fabs(y - f_surf(x,y,0,0)) < 0.5*diam && level < maxlevel+1 );
		//refine( fabs(x) < 0.4*L0  && fabs(y) < L0/pow(2, 0.5*(minlevel+maxlevel)) && level < maxlevel );
		//refine( fabs(x) > 0.4*L0  && fabs(y-(-5.*fabs(x)+2.0*L0)) < L0/pow(2, 0.5*(minlevel+maxlevel)) && level < maxlevel );
		substrate (fs);
		jet_in(f, fs);
	}
	f.refine = f.prolongation = fraction_refine;
	fs.refine = fs.prolongation = fraction_refine;
	fprintf(ferr, "v_in=%g", v_in);
	event("vtk_file");

	/**
	We initialize the reference velocity. */
}

/**
The gravity vector is aligned with the channel and viscosity is
unity. */
#if !REDUCED
event acceleration (i++) {
	face vector av = a;
	foreach_face(y)	av.y[] += grav*f[];//m^2/s
}
#endif
/**
We check for a stationary solution. */

event snapshot (t += 0.01) {
	char name[80];
	sprintf (name, "snapshot-%g", t);
	scalar pid[];
	foreach()	pid[] = fmod(pid()*(npe() + 37), npe());
			boundary ({pid});
	dump (name);
}
scalar un[];
#undef SEPS
#define SEPS 1e-15
event logfile (i+=10)
{
	double avg = normf(u.x).avg, du = change (u.x, un)/(avg + SEPS);
	fprintf (ferr, "%d %d %g %g %d %d %d %d %d %d %.3g %.3g %.3g %.3g %.3g\n",
	maxlevel, i, t, dt,
	mgp.i, mgp.nrelax, mgp.minlevel,
	mgu.i, mgu.nrelax, mgu.minlevel,
	du, mgp.resa*dt, mgu.resa, statsf(u.x).sum, normf(p).max);
}
//event pressure_correction(i++){
//	norm normp = normf (p);
//	norm normpf = normf (pf);
//	fprintf(ferr, "pressure_avg=%g pf=%g\n", normp.avg, normpf.avg);
//	foreach() {p[] -= normp.avg; pf[] -=normpf.avg;}
//}
//Output
//event vtk_file (i+=10; t<100){
event vtk_file (t+=0.01; t<2.25){
	char subname[80]; sprintf(subname, "hjump");
	scalar l[];
	vorticity (u, omega);
	foreach() {l[] = level; omega[] *= 1 - fs[];}
	output_vtu_MPI( (scalar *) {fs, f, omega, p, l}, (vector *) {u, uf, a}, subname, 0.0);
}

event adapt (i++) {
	adapt_wavelet ({f, fs, u}, (double[]){1e-3, 1e-4, 1e-2, 1e-2, 1e-2}, maxlevel=maxlevel, minlevel=minlevel);
	fs.refine = fs.prolongation = fraction_refine;
	boundary({fs});
}
