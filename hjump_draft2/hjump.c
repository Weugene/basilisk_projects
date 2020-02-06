/**
# Stokes flow through a complex porous medium


This tests mainly the robustness of the representation of embedded
boundaries and the convergence of the viscous and Poisson
solvers. */
#define BRINKMAN_PENALIZATION 1
#define DEBUG_BRINKMAN_PENALIZATION 1
#define REDUCED 1
#include "../src_local/centered-weugene.h"
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
#include "tension.h"
#if REDUCED
#include "reduced.h"
#endif
#include "view.h"
#include "../src_local/output_vtu_foreach.h"
#define FILTERED
/**
We will vary the maximum level of refinement, starting from 5. */

int maxlevel = 10;
int minlevel = 4;
#define diam 1.2e-3//meter 5e-3
#define Q 3.33e-5//cm^3/s Q=v*S
#define v_in 4.0*Q/(pi*sq(diam))//inlet velocity
double hjet = 3e-2;//meter
double length = 0.5e-2;//4.1e-3
double Tac = 1e-1;

double SIGMA = 30e-3;//30-70e-3
double MUL = 3.5e-3, MUG = 18.5e-5;
double RHOL = 1e+3, RHOG = 1.2;
scalar f0[], fs[], omega[];

//u.n[top]  = f0[] > 0 ? dirichlet( - v_in*f0[]*min(t/Tac, 1)) : neumann(0);
//u.t[top]  = f0[] > 0 ? dirichlet(0) : neumann(0);
//uf.n[top]  = f0[] > 0 ? dirichlet( - v_in*f0[]*min(t/Tac, 1)) : neumann(0);
//uf.t[top]  = f0[] > 0 ? dirichlet(0) : neumann(0);
u.n[top]  = dirichlet( - v_in*f0[]*min(t/Tac, 1));
u.t[top]  = dirichlet(0);
uf.n[top]  = dirichlet( - v_in*f0[]*min(t/Tac, 1));
uf.t[top]  = dirichlet(0);
p[top]    = neumann(0);
pf[top]    = neumann(0);
f[top]    = f0[] > 0 ? dirichlet( f0[]) : neumann(0);

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
p[left]   = neumann(0);
pf[left]   = neumann(0);
f[left]    = neumann(0);

u.n[right] = neumann(0);
u.t[right] = neumann(0);
uf.n[right] = neumann(0);
uf.t[right] = neumann(0);
p[right]   = neumann(0);
pf[right]   = neumann(0);
f[right]    = neumann(0);
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
//	return 0.5*(L0 - 2.0*hjet)*(-1.0 - tanh(100.0*(sq(x) - sq(0.37*L0))/sq(L0)));
//	return pow(pow(R,n) -pow(x,n), 1./n) - R;
	return ((fabs(x)<0.4*L0 && y<=0) ? 1 : (y<=5.*fabs(x)-2.0*L0) ? 1 : -1);
}
void substrate (scalar fs)
{
//	face vector ffs[];
	vertex scalar phi[];
	foreach_vertex()
	phi[] = ((fabs(x)<0.4*L0 && y<=0) ? 1 : (fabs(x)>0.4*L0 && y<=-5.*fabs(x)+2.0*L0) ? 1 : 0);
//	boundary ({phi});
	foreach() fs[] = 0.25*(phi[] + phi[1,0] + phi[0,1] + phi[1,1]);
//	fraction (fs, (f_surf(x, y, 0.26, 10)));
//	fraction (fs, (  4*     -0.5*(L0 - 2.0*hjet)*(1.0 + tanh(50.0*(x*x - 0.37*L0*0.37*L0)/(L0*L0))) - y));
//	foreach() fprintf(ferr, "fs=%g", fs[]);
//	fraction (fs, y - f_surf(x));
}

double f_jet(double x, int n){
	return hjet - length*pow(1 -pow(x/(0.5*diam),n), 1./n);
}

void jet_in (scalar f1,scalar fs1)
{
	scalar fc[];
	double R = (sq(diam) + sq(0.4*L0))/(2.0*diam);
	double xc = 0, yc = -R + diam;
//	fraction (f0, sq(0.5*diam) - sq(x));//sq(0.5*diam) - sq(y - hjet) - sq(x)
//	foreach() f1[] = f0[]*(y > hjet-length);//0.5*hjet
	fraction (f0, 1 - pow((y-hjet)/length,10) - pow(x/(0.5*diam),10) );//sq(0.5*diam) - sq(y - hjet) - sq(x)
	fraction (fc, (sq(R) - sq(x - xc) - sq(y - yc) > 0 && y > 0)? 1 : -1 );
	foreach() {f0[] =clamp(f0[] + fc[], 0, 1); f1[] = f0[];} //clamp(f0[]+fs1[], 0, 1);
}

int main()
{
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
	G.y = -9.8;
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
				u.y[] = - v_in*clamp(f[]-fs[], 0, 1)*min(t/Tac, 1);
			}
			boundary (all); // this is necessary since BCs depend on embedded fractions
		}while (adapt_wavelet({fs,f}, (double []){1e-6, 1e-4}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
		refine( fabs(x) < 0.4*L0  && fabs(y) < L0/pow(2, 0.5*(minlevel+maxlevel)) && level < maxlevel );
		refine( fabs(x) > 0.4*L0  && fabs(y-(-5.*fabs(x)+2.0*L0)) < L0/pow(2, 0.5*(minlevel+maxlevel)) && level < maxlevel );
		substrate (fs);
		jet_in(f, fs);
	}
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
	foreach_face(y)	av.y[] += -9.8*f[];//m^2/s
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
//event vtk_file (i+=1; t<100){
event vtk_file (t+=0.001; t<100){
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

//	refine( fabs(x) < 0.4*L0  && fabs(y) < L0/pow(2, 0.5*(minlevel+maxlevel)) && level < maxlevel );
}
