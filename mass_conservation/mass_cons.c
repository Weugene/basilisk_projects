#define BRINKMAN_PENALIZATION 1
#define DEBUG_BRINKMAN_PENALIZATION 1
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

int maxlevel = 8;
int minlevel = 5;
int Nobst = 2; //250
scalar fs[], omega[];
double U0=0.01, rhol=1e+3, sig=73e-3, Lchar=5e-3, mul=1e-3, grav=9.8;
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
#if BRINKMAN_PENALIZATION == 4
    lambda_slip = L0*pow(2.0, -maxlevel-1);
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

void calc_norms(scalar fs, vector n_sol, int ns){
    double xc[ns], yc[ns], R[ns];
    double dist = L0/ns, magn2;
    double size = 0.25*dist;
    for (int i = 0; i < ns; i++) {
        xc[i] = 0.25*L0;
        yc[i] = -0.5*L0 + 0.5*dist + dist*i;
        R[i] = size;
        fprintf(ferr, "xcycR %g %g %g \n", xc[i], yc[i], R[i]);
    }
    trash ({n_sol});
    foreach() {
        if (fs[]>0){
            for(int i=0; i<ns; i++) {
                magn2 = sq(x - xc[i]) + sq(y - yc[i]) + 1e-15;
                if (magn2 <= sq(R[i])) {
                    n_sol.x[] = (x - xc[i]) / sqrt(magn2);
                    n_sol.y[] = (y - yc[i]) / sqrt(magn2);
                    fprintf(ferr, "n %g %g %g %g %g \n", x, y, n_sol.x[], n_sol.y[], magn2);
                }
            }
        }else{
            foreach_dimension() n_sol.x[] = 0;
        }
    }
}

#if BRINKMAN_PENALIZATION==4
event properties(i++){
    calc_norms(fs, n_sol, Nobst);
//    foreach() {
//        if (fs[] > SEPS && fs[] < 1 - SEPS) {
//            n_sol.x[] = (fs[] - fs[-1]) / Delta;
//            n_sol.y[] = (fs[] - fs[0, -1]) / Delta;
//            mag_n = sqrt(sq(n_sol.x[]) + sq(n_sol.y[]));
//            n_sol.x[] /= (mag_n + SEPS);
//            n_sol.y[] /= (mag_n + SEPS);
//        }else{
//            n_sol.x[] = 0.0;
//            n_sol.y[] = 0.0;
//        }
//    }
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

event logfile (t += 0.01)
{
	double avggas = sq(L0) - normf(f).avg;
	scalar umag[];
	foreach() umag[] = norm(u); // the length of u
	norm statu = normf_weugene(umag, fs); // outputs avg, rms, max, volume
	fprintf (ferr, "%d %d %g %g %g %g %g %g \n", maxlevel, i, t, dt, avggas, statu.avg, statu.rms, statu.max);
}

void correct_press(scalar p, int i){
    double press = 0;
    int ip = 0;
//    FILE *fp1;
//    char subname[80]; sprintf(subname, "nameYouWant-%d", pid());
//    fp1 = fopen(subname, "a");
#if 1 // Left bottom Corner
    foreach(){
        if (ip == 0){
            press = p[];
            ip++;
            @if _MPI
                MPI_Bcast(&press, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            @endif
        }
    }
#else //average value
    press = normf(p).avg;
#endif
    foreach(){
        p[] -= press;
    }
//    fprintf(ferr, "p %g \n", press);
//    fclose(fp1);
}

event end_timestep(i++){
    correct_press(p, i);
}
//Output
//event vtk_file (i++; t<10){
event vtk_file (t += 0.01; t<10){
	char subname[80]; sprintf(subname, "mc");
	scalar l[], npid[];
	vorticity (u, omega);
//	vector n_ss[];
//    calc_norms(fs, n_ss, Nobst);
	foreach() {l[] = level; omega[] *= 1 - fs[]; npid[] = pid();}
    scalar divu[];
    foreach() {
        divu[]=0;
        foreach_dimension() divu[] += (uf.x[1]-uf.x[])/Delta;
    }
//#ifndef DEBUG_BRINKMAN_PENALIZATION && BRINKMAN_PENALIZATION == 4
	output_vtu_MPI( subname, (iter_fp) ? t + dt : 0, (scalar *) {fs, f, omega, p, l, npid, divu}, (vector *) {u, a});
//#else
//   output_vtu_MPI( (scalar *) {fs, f, omega, p, l, npid}, (vector *) {u, a, n_sol, target_U, dbp, total_rhs, utau, grad_utau_n}, subname, 1 );
//#endif
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
