#define BRINKMAN_PENALIZATION 1
//#define DEBUG_BRINKMAN_PENALIZATION 1
#define DEBUG_MINMAXVALUES 1
#define FILTERED
#define JACOBI 1
#undef SEPS
#define SEPS 1e-30
scalar f0[], divu[];
face vector fs_face[];

#include "grid/octree.h"
#include "../src_local/centered-weugene.h"
//#include "two-phase.h"
#include "../src_local/three-phase-weugene.h"
//#include "tension.h"
//#include "navier-stokes/conserving.h"
/**
We also need to compute distance functions (to describe the solid geometry
), use visualisation functions. */

#include "distance.h"
#include "view.h"

/**
On supercomputers we need to control the maximum runtime and we check
performances. */

#include "maxruntime.h"
#include "navier-stokes/perfs.h"
#include "../src_local/output_vtu_foreach.h"
/**
## Importing the geometry 

This function computes the solid fraction given a pointer to an STL
file, a tolerance (maximum relative error on distance) and a 
maximum level. */

void fraction_from_stl (scalar f, FILE * fp, double eps, int maxlevel){
	/**
	We read the STL file and compute the bounding box of the model. */
	coord * p = input_stl (fp);
	coord min, max;
	bounding_box (p, &min, &max);
	fprintf(ferr, "min= (%g %g %g), max = (%g %g %g) \n", min.x, min.y, min.z, max.x, max.y, max.z);
	double maxl = -HUGE;
	foreach_dimension() if (max.x - min.x > maxl)
	maxl = max.x - min.x;
  
	/**
	We initialize the distance field on the coarse initial mesh and
	refine it adaptively until the threshold error (on distance) is
	reached. */

	scalar d[];
	distance (d, p);
	while (adapt_wavelet ({d}, (double[]){eps*maxl}, maxlevel, 5).nf);

	/**
	We also compute the volume fraction from the distance field. We
	first construct a vertex field interpolated from the centered field
	and then call the appropriate VOF functions. */

	vertex scalar phi[];
	foreach_vertex()
		phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1]
				+ d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.; //left bottom vertex
	fractions (phi, f);
}

int minlevel = 5, maxlevel = 9;
int level = 6;
double Ldomain = 1.02, Lch = 1.02;//6e-3
double uemax = 0.1;
double rhol = 1, rhog = 0.1, mul = 1, mug = 1e-1;
double sig = 0.0005, pdelta = 1;
//double rhol = 1600, rhog = 10, mul = 1, mug = 1e-3;
//double sig = 0.074, pdelta = 330.0;
double U0, RE, CA, LA, WE;
/**
We need additional (fraction) fields for the composite geometry and for the
(inflow) boundary condition. */



int main (int argc, char * argv[]) {

	maxruntime (&argc, argv);
	if (argc > 1) maxlevel = atoi(argv[1]); //convert from string to int
	fprintf(ferr, "maxlevel = %d \n", maxlevel);
	stokes = true;
	//periodic(top);
	//periodic(front);
	init_grid (1 << level);
    	U0 = pdelta*sq(Lch/10.0)/(2.0*mul*Lch);
    	RE = rhol*U0*Lch/mul;
    	CA = mul*U0/sig;
    	LA = sig*rhol*Lch/sq(mul);
    	WE = rhol*sq(U0)*Lch/sig;
	rho1 = 1.0; // water
	rho2 = rhog/rhol; // air
	rho3 = max(rho1, rho2); // air
	mu1 = 1.0/RE;
	mu2 = mug/mul/RE;
	mu2 = max(mu1, mu2);
//	f.sigma = 1.0/WE;
    	fprintf(ferr, "Tension.h module is switched off");
    	eta_s= 1e-6;
	TOLERANCE = 1e-3;
	NITERMAX = 15;
	size (Ldomain);
    	double sh=0.0*L0;
	origin (-sh,-sh,-sh);
	fprintf(ferr, "U0=%g dp=%g Lch=%g sigma=%g"
                  "mu1=%g mu2=%g rho1=%g rho2=%g"
                  "Re=%g Ca=%g We=%g La=%g", U0, pdelta, Lch, sig,
                   mu1, mu2, rho1, rho2, RE, CA, WE, LA);
	/**
	We need to tell the code that both `fs` and `f0` are volume
	fraction fields. */
	for (scalar s in {fs,f0}) s.refine = s.prolongation = fraction_refine;

	run();
}

/**
## Boundary conditions
The inflow condition */
u.n[left] = neumann(0);//dirichlet(U0*(1-fs[]));
u.t[left] = neumann(0);
u.r[left] = neumann(0);
uf.n[left] = neumann(0);
uf.t[left] = neumann(0);
uf.r[left] = neumann(0);
p[left]   = dirichlet(pdelta);
pf[left]  = dirichlet(pdelta);
f[left]   = dirichlet(1.0 - fs[]);
/*The outflow condition */
u.n[right] = neumann(0);
u.t[right] = neumann(0);
u.r[right] = neumann(0);
uf.n[right] = neumann(0);
uf.t[right] = neumann(0);
uf.r[right] = neumann(0);
p[right]   = dirichlet(0);
pf[right]  = dirichlet(0);
f[right]   = neumann(0);


/**
## Initial conditions

We can optionally restart, otherwise we open the STL file and
initialize the corresponding fraction. We also initialize the `f0`
field used for the inflow condition and set the initial epoxy level
and velocity field. */

event init (t = 0) {
	if (!restore (file = "restart")) {
		if (pid()>0) {
			fprintf(ferr, "you have run in parallel reading ONLY in serial N parallel=%d.. "
				 "Run in serial save dump, then rename it to restart, then run in parallel", pid());
			exit(1992);
		}
        	FILE * fp = fopen ("cube.stl", "r");
		fprintf(ferr, "opened. \n");
		fraction_from_stl (fs, fp, 1e-3, maxlevel);
		fprintf(ferr, "stl saved in fs \n");
//		face_fraction (fs, fs_face);
        	foreach_face() fs_face.x[] = 0.5*(fs[-1] + fs[]);
        	boundary((scalar *){fs_face});

		fclose (fp);
//        	foreach() {f[] += fs[]; f0[] = f[];}
		refine(x < X0 + Ldomain/pow(2., minlevel) && level < maxlevel);
		boundary ({fs, f, f0, u});
		DT = 1e-9;
		event("snapshot");
	}
	??DT = 1e-6;
}

event set_dtmax (i++) if (i<500) DT *= 1.05;



event end_timestep (i+=100){
	double avggas = sq(L0) - normf(f).avg;
	foreach() {
		divu[] = 0;
		foreach_dimension() divu[] += (uf.x[1] - uf.x[])/Delta;
	}
	double Linf_u = -10;
	foreach( reduction(max:Linf_u) ){
		if (fabs(divu[]) > Linf_u) Linf_u = fabs(divu[]);
	}
	if (pid()==0) fprintf (ferr, "i=%d t=%g dt=%g iter_p=%d iter_u=%d AvgGas=%g divu=%g \n", i, t, dt, mgp.i, mgu.i, avggas, Linf_u);
	if (pid()==0) fprintf (stderr, "%d %g %g %d %d\n", i, t, dt, mgp.i, mgu.i);
	//double eps_arr[] = {1,1,1, 1,1,1, 1,1,1, 1,1,1};
    //MinMaxValues({f,fs,p,u,uf,dbp}, eps_arr);
	//double eps_arr[] = {1,1,1};
    //MinMaxValues({f,fs,p}, eps_arr);
}
/**
## Animations

We generate animations of the composite surface (as represented by the
solid fraction fs) and of the gas-liquid interface, colored with the
height.
The computations above were done on the Irene supercomputer using 12
levels of refinement. */
event movie (t += 0.01; t <= 10) {
//event movie (i += 100) {
    view (fov = 50, camera="iso",
        tx = 0, ty = 0.2,
        width = 1024, height = 768);
	clear();
	draw_vof ("fs", fc = {0.5,0.5,0.5});
	//scalar umag[];
	//foreach() umag[] = norm(u);
	draw_vof ("f", fc = {0,0,1});
	save ("movie.mp4");

//	draw_vof ("fs", fc = {0.5,0.5,0.5});//fs is grey
//	draw_vof ("f", color = "Z", min = -0.1, max = 0.1, linear = true);
//
//	lambda2 (u, l2);
//	isosurface ("l2", -100);
//	save ("l2.mp4");
}

#if DUMP
event snapshot (i +=100) {
	char name[80];
	sprintf (name, "dump-%d", i);
	dump (file = name);
	fprintf(ferr,"dumped file=%s in t=%g dt=%g i=%d \n", name, t, dt, i);
}
#endif

//Output
//event vtk_file (t += 0.01){
event vtk_file (i += 1){
    char subname[80]; sprintf(subname, "comp");
    scalar l[];
    //vorticity (u, omega);
    foreach() {l[] = level;}
    output_vtu_MPI( (scalar *) {f, fs, l}, (vector *) {u}, subname, 0);
}

/**
## Mesh adaptation

This computation is only feasible thanks to mesh adaptation, based
both on volume fraction and velocity accuracy. */
#define ADAPT_SCALARS {f, fs, u}
#define ADAPT_EPS_SCALARS {1e-4,1e-4,uemax,uemax,uemax}
event adapt (i++){
    double eps_arr[] = ADAPT_EPS_SCALARS;
    MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
	fs.refine = fs.prolongation = fraction_refine;
	boundary({fs});
}

event stop(t = 10);
