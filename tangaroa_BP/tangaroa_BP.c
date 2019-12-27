/**
# Two-phase flow around RV Tangaroa

This is an improved version of the [famous Gerris
example](http://gerris.dalembert.upmc.fr/gerris/examples/examples/fs.html),
illustrating the combination of complex solid boundaries, air-water
turbulent flows and reduced gravity approach.

We use the centered Navier--Stokes solver, two-phase flow and the
momentum-conserving option. Note that the momentum-conserving option
is crucial to obtain stable solutions for this air-water density ratio
configuration. */

#define BRINKMAN_PENALIZATION 1
#define DEBUG_BRINKMAN_PENALIZATION 1
#define DEBUG_MINMAXVALUES 1

#undef SEPS
#define SEPS 1e-30


#include "grid/octree.h"
#include "../src_local/centered-weugene.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"

/**
We also need to compute distance functions (to describe the ship
geometry), use reduced gravity and visualisation functions. */

#include "distance.h"
#include "reduced.h"
#include "view.h"
#include "lambda2.h"

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

/**
## Main function

We can change both the maximum level of refinement and the [Froude
number](https://en.wikipedia.org/wiki/Froude_number) at runtime.

[RV Tangaroa](https://en.wikipedia.org/wiki/RV_Tangaroa) is 70 metres
long. If we assume that it moves at 20 knots (twice its actual cruise
speed), this gives a Froude number of approx 0.4. */

int MAXLEVEL = 12;
int LEVEL = 9;
double FROUDE = 0.4;
double uemax = 0.1; 
/**
We need additional (fraction) fields for the ship geometry and for the
(inflow) boundary condition. */

scalar fs[], f0[], l2[], omega[];

int main (int argc, char * argv[]) {

	maxruntime (&argc, argv);
	if (argc > 1) MAXLEVEL = atoi(argv[1]); //convert from string to int
	fprintf(ferr, "maxlevel = %d", MAXLEVEL);
	if (argc > 2) FROUDE = atof(argv[2]);
	fprintf(ferr, "Froude = %g", FROUDE);
	init_grid (32);
	rho1 = 1.; // water
	rho2 = 1./815.; // air
    eta_s= 1e-15;
	/**
	The length of the ship is unity and the domain is five times
	larger. We change the origin so that the ship is not too close to
	the inflow. */

	size (5.);
	origin (-L0/2.,-L0/3.,-L0/2.);

	/**
	We need to tell the code that both `fs` and `f0` are volume
	fraction fields. */

	for (scalar s in {fs,f0}) {
		s.refine = s.prolongation = fraction_refine;
	}

	/**
	Since the ship length is one and the velocity one, the acceleration
	of gravity is...*/
	G.z = - 1./sq(FROUDE);
	run();
}

/**
## Boundary conditions

The inflow condition fixes the velocity (unity) and the water level
(using `f0`). */
u.n[bottom] = dirichlet(1);
p[bottom]   = neumann(0.);
pf[bottom]  = neumann(0.);
f[bottom]   = f0[];

/**
Outflow uses standard Neumann/Dirichlet conditions.  */
u.n[top]  = neumann(0.);
p[top]    = dirichlet(0.);
pf[top]   = dirichlet(0.);

/**
Boundary conditions for the solid and fraction tracers. */

fs[back] = 0;
f[back]  = 1;

/**
Not sure whether this is really useful. */
uf.n[left] = 0.;
uf.n[right] = 0.;
/**
## Initial conditions

We can optionally restart, otherwise we open the STL file and
initialize the corresponding fraction. We also initialize the `f0`
field used for the inflow condition and set the initial water level
and velocity field. */

event init (t = 0) {
	if (!restore (file = "restart")) {
		FILE * fp = fopen ("tangaroa.stl", "r");
		fprintf(ferr, "open");
		fraction_from_stl (fs, fp, 5e-4, MAXLEVEL);
		fprintf(ferr, "stl");
		fclose (fp);
		fraction (f0, - z);
		fprintf(ferr, "fraction");
		foreach() {
			f[] = f0[];
			u.y[] = 1.;
		}
		fprintf(ferr, "foreach");
		boundary ({f,u.y});
	}
}

event logfile (i++){
	if (pid()==0) fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
	//double eps_arr[] = {1,1,1, 1,1,1, 1,1,1, 1,1,1};
    //MinMaxValues({f,fs,p,u,uf,dbp}, eps_arr);
	//double eps_arr[] = {1,1,1};
    //MinMaxValues({f,fs,p}, eps_arr);

}
/**
## Animations

We generate animations of the ship surface (as represented by the
solid fraction) and of the air-water interface, colored with the
height.

Several classical features of ship wakes are recognisable: breaking
bow wave, breaking stern divergent wave, turbulent boundary layer,
Kelvin waves etc...

![Evolution of the air-water interface](fs/movie.mp4)(width="800" height="600")

We also use the $\lambda_2$ criterion to display the turbulent
vortical structures generated by the airflow. The air recirculation at
the top of the steep primary Kelvin waves is particularly noticeable.

![Turbulent vortical structures](fs/l2.mp4)(width="800" height="600")

The computations above were done on the Irene supercomputer using 12
levels of refinement. */
event movie (t += 100.01; t >= 10) {
    view (fov = 5.86528,
        quat = {0.515965,0.140691,0.245247,0.808605},
        tx = -0.07438, ty = -0.0612925,
        width = 1024, height = 768);
	clear();
	draw_vof ("fs");
	scalar Z[];
	Z[back] = dirichlet (z);
	foreach() Z[] = z;
	boundary ({Z});
	draw_vof ("f", color = "Z", min = -0.1, max = 0.1, linear = true);
	save ("movie.mp4");

	draw_vof ("fs", fc = {0.5,0.5,0.5});//fs is grey
	draw_vof ("f", color = "Z", min = -0.1, max = 0.1, linear = true);

	lambda2 (u, l2);
	isosurface ("l2", -100);
	save ("l2.mp4");
}

#if DUMP
event snapshot (t +=0.01) {
	char name[80];
	sprintf (name, "dump-%d", i);
	lambda2 (u, l2);
	dump (file = name);
}
#endif

//Output
//event vtk_file (t += 0.01){
//    char subname[80]; sprintf(subname, "br");
//    scalar l[];
//    vorticity (u, omega);
//    foreach() {l[] = level; omega[] *= 1 - fs[];}
//    output_vtu_MPI( (scalar *) {l, omega, fs, p, l2}, (vector *) {u, uf, dbp}, subname, 0);
//}

/**
## Mesh adaptation

This computation is only feasible thanks to mesh adaptation, based
both on volume fraction and velocity accuracy. */
#define ADAPT_SCALARS {f, fs, u}
#define ADAPT_EPS_SCALARS {0.01,0.01,uemax,uemax,uemax}
event adapt (i++){
    double eps_arr[] = ADAPT_EPS_SCALARS;
//    MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = MAXLEVEL, minlevel = 5);
//    calc_solid(fs, n_sol, target_U);
}
