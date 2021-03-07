#define BRINKMAN_PENALIZATION 1
#define DEBUG_BRINKMAN_PENALIZATION 1
#define DEBUG_OUTPUT_VTU_MPI
#define FILTERED
#include "../src_local/centered-weugene.h"
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
#include "tension.h"
#include "view.h"
#include "../src_local/output_vtu_foreach.h"

int maxlevel = 10;
int minlevel = 2;
scalar fs[], omega[], divu[];
double U0, Umax, rhol = 1, sig, mul = 1, deltap;
double RE=0.5, CA=0.01, FR, B;
double Rrho = 1, Rmu = 100;
double Hch = 1, Ldomain = 8, Radius_b = 0.3;
double deltat=0.05, timef=100;
#define CASE 2

#if CASE==1
double Zet = 7.4;
#else
double Radius_b2 = 1.2;
#endif
#define uy (deltap/(2*mu1*Ldomain))*(sq(0.5*Hch) - sq(y))*(1-fs[])
u.n[left]  = dirichlet(uy);
u.t[left]  = dirichlet(0.);
uf.n[left] = dirichlet(uy);
uf.t[left] = dirichlet(0.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(1);

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
uf.n[right] = neumann(0.);
uf.t[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);
f[right]   = dirichlet(1);

void obstacles (scalar fs)
{
	fraction (fs, fabs(y) - 0.5*Hch  );//
}

#if CASE==1
	#define func_bubble  (sq(x) + sq(y) - sq(Radius_b))
#else
	#define func_bubble (pow(fabs(x)/Radius_b2, 10) + pow(fabs(y)/Radius_b, 10) - sq(1))
#endif
void bubbles (scalar f)
{
	#if CASE==1
    	fraction (f, sq(x) + sq(y) - sq(Radius_b) );
    #else
    	fraction (f, pow(fabs(x)/Radius_b2, 10) + pow(fabs(y)/Radius_b, 10) - sq(1) );
    #endif

	vertex scalar phi[];
	face vector ff[];
	foreach_vertex() {
		phi[] = HUGE;
		phi[] = intersection (phi[], func_bubble);
	}
	boundary ({phi});
	fractions (phi, f, ff);
}


int main(int argc, char * argv[])
{
    if (argc > 1) {
        Radius_b = atof(argv[1]); 
    }
    #if CASE==1
    if (argc > 2) {
        Zet = atof(argv[2]); 
    }
    #else
    if (argc > 2) {
        CA = atof(argv[2]); 
    }
    #endif
	if (argc > 3) {
		RE = atof(argv[3]);
	}
    if (argc > 4) {
        maxlevel = atoi(argv[4]);
    }
    if (argc > 5) {
        Rmu = atof(argv[5]);
    }
	eta_s = 1e-6;
	TOLERANCE = 1e-4;
	NITERMAX = 30;
	N = 1 << 6;

	rho1 = rhol; rho2 = rho1/Rrho;
	mu1 = mul; mu2 = mu1/Rmu;
	#if CASE==1
    Hch = 1.0;
    U0 = 1.0;
    deltap = 12.0*mu1*Ldomain*U0/sq(Hch);
    B = deltap/(2.0*mu1*Ldomain);
    sig = Zet*sq(Radius_b)*B*mu1;
    RE = rho1*U0*Hch/(2.0*mu1);
    CA = mu1*1.5*U0/sig;
    fprintf(ferr, "Zet=%g ", Zet);
    #else
    Hch = 1.0;
    Radius_b = 0.45*Hch;
    Radius_b2 = 2.0*Radius_b; //ratio of drop length to radius = 4
    U0 = 2.0*RE*mu1/(rho1*Hch);
    deltap = 12.0*mu1*Ldomain*U0/sq(Hch);
    B = deltap/(2.0*mu1*Ldomain);
	sig = mu1*1.5*U0/CA;
	#endif
    timef = Ldomain/U0;
    deltat = timef/500.0;
    size (Ldomain);
	origin (-6*Radius_b, -0.5*L0);
	f.sigma = sig;
	fprintf(ferr, "RE=%g CA=%g \n"
			   "deltap=%g U0=%g B=%g \n"
			   "mu1=%g mu2=%g rho1=%g rho2=%g \n"
	           "sigma=%g Rb=%g L0=%g\n",
	           RE, CA, deltap, U0, B, mu1, mu2, rho1, rho2, f.sigma, Radius_b, L0);
	run();
}

scalar un[];

event init (t = 0) {
	if (!restore (file = "restart")) {
		int it = 0;
		do {
			obstacles (fs);
			bubbles(f);
			foreach() u.x[] = uy;
		}while (adapt_wavelet({fs, f, u}, (double []){1e-3, 1e-3, 1e-3, 1e-3}, maxlevel=maxlevel, minlevel=2).nf != 0 && ++it <= 10);
		DT=1e-10;
	}
	event("vtk_file");
}

event set_dtmax (i++) if (i<500) DT *= 1.05;

//event snapshot (t += 0.5; t <= 10.8) {
//	char name[80];
//	sprintf (name, "snapshot-%g", t);
//	scalar pid[];
//	foreach()	pid[] = fmod(pid()*(npe() + 37), npe());
//	boundary ({pid});
//	dump (name);
//}

event logfile (i +=100)
{
	double avggas = L0*Hch - normf(f).avg;
	scalar umag[];
	double xcg = 0, volume = 0, volumeg = 0, velamean = 0, velgx = 0, velgy = 0, dvtmp;
	foreach(reduction(+:xcg) reduction(+:volume) reduction(+:volumeg)
			reduction(+:velgx) reduction(+:velgy)
			reduction(+:velamean) ) {
		if (fs[]<1){
			dvtmp = (1.0 - f[])*(1.0 - fs[])*dv();
			volumeg += dvtmp;//gas liquid
			volume += (1.0 - fs[])*dv();//gas liquid
			umag[] = norm(u); // the length of u
			xcg   += x*dvtmp;// Along x
			velamean += (1.0 - fs[])*umag[]*dv();//mean velocity of gas
			velgx += u.x[]*dvtmp;//mean velocity of gas Ox
			velgy += u.y[]*dvtmp;//mean velocity of gas Oy
		}
	}
	xcg /= volumeg; velgx /= volumeg; velgy /= volumeg; velamean /= volume;
	//norm statu = normf_weugene(umag, fs); // outputs avg, rms, max, volume
	fprintf (ferr, "%d %d %d %g %g %g %g %g %g %g %g %g %g %g \n", maxlevel, i, iter_fp, t, dt, avggas, velgx, velgy, velamean, velgx/(sq(Radius_b)*B), velgx/U0, (velgx/U0 - 1), log(fabs(velgx/U0 - 1)),  xcg);
    if (xcg > Ldomain-Radius_b2) {
        fprintf(ferr, "a drop is outside of the domain");
        exit(197);
    }
}
event time_end(i++){
	double divu_min = 1e10, divu_max = -1e10;
	foreach(reduction(min:divu_min) reduction(max:divu_max)) {
		divu[]=0;
		foreach_dimension() divu[] += (uf.x[1]-uf.x[])/Delta;
		if(divu[] > divu_max) divu_max = divu[];
		if(divu[] < divu_min) divu_min = divu[];
	}
	fprintf (ferr,"divu_min=%g divu_max=%g \n", divu_min, divu_max);
}
//Output
event vtk_file (t += deltat; t < timef){
	char subname[80]; sprintf(subname, "bubble_in_channel");
	scalar l[];
	vorticity (u, omega);
	foreach() {l[] = level; omega[] *= 1 - fs[]; }

	output_vtu_MPI( (scalar *) {fs, f, omega, p, l, divu}, (vector *) {u, a}, subname, 0 );
}

#define ADAPT_SCALARS {f, fs, u.x, u.y}
#define ADAPT_EPS_SCALARS {1e-3, 1e-3, 0.01*U0, 0.01*U0}
event adapt (i++){
	double eps_arr[] = ADAPT_EPS_SCALARS;
//	MinMaxValues(ADAPT_SCALARS, eps_arr);
	adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
	for (scalar s in {f,fs})
        s.refine = s.prolongation = fraction_refine;
	boundary({f,fs});	
}



#if 0
#!/usr/bin/gnuplot -c
print "args:".ARG1
if (strlen(ARG1) == 0) print "Usage: " . ARG0 . " picture case. 1-relative Err. 2-average velocity, 3-max velocity"; exit
# Output W3C Scalable Vector Graphics
##set terminal pdf
set terminal postscript eps enhanced color font 'Helvetica,10'
set grid
set key spacing 2
set key top left
set tics font "Helvetica,10"
set for [i=1:7] linetype i dt i
set style line 1 lt 1 lc rgb "blue" lw 3 pt 1 ps 2
set style line 2 lt 1 lc rgb "red" lw 3 pt 2 ps 2
set style line 3 lt 1 lc rgb "green" lw 3 pt 4 ps 2
set style line 4 lt 1 lc rgb "green" lw 3 pt 33 ps 1
set style line 5 lt 1 lc rgb "blue" lw 3 pt 10 ps 1
set style line 6 lt 1 lc rgb "red" lw 3 pt 12 ps 1
set style line 7 lt 3 lc rgb "black" lw 3 pt 2 ps 1
mcase=ARG1+0
array A[3]
array N[3]
array Title[3]
array MAXFi[3]
array lstyles[3]
A[1]=0.08
A[2]=0.125
A[3]=0.2
tmax=0.4
dataname(n) = sprintf("m%g",n)
set xlabel 'time'  font ",10"

do for [i=1:3] {
    N[i]=sprintf("m%g",A[i]);
    Title[i]=sprintf("r1=%g",A[i]);
    stats  N[i] u 3:5 name "XX";
    MAXFi[i]=XX_max_y;
    lstyles[i]=i;
}
set xr [0:tmax]

if (mcase==1){
    set yr [-1e-5:1e-5];
    set format y "10^{%L}"
    set ylabel 'Relative Error, (fi-fi0)/fi0'  font ",10"
    plot for[i=1:3] N[i] u 3:(($5 - MAXFi[i])/MAXFi[i]) t Title[i] w lp ls i;
}

if (mcase==2){
    set ylabel 'Average velocity'  font ",10"
    plot for[i=1:3] N[i] u 3:6 t Title[i]  w lp ls i
}

if (mcase==3){
    set ylabel 'Max velocity'  font ",10"
    plot for[i=1:3] N[i] u 3:8 t Title[i] w lp ls i
}
#endif
