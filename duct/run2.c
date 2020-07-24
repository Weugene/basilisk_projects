scalar my_residual[], divutmp[];
#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "../src_local/utils-weugene.h"
#include "../src_local/output_vtu_foreach.h"
#include "view.h"
#include "maxruntime.h"

//Channel cross section Lyy*Lzz
#define Lyy 1.
#define Lzz 1.
#define Uin 1. //U inlet
#define Re 2.5 //Reynolds

int LEVEL = 8;
scalar cs[], un[];
face vector muv[];

int main (int argc, char * argv[]) {
  maxruntime (&argc, argv);
  if (argc > 1)
    LEVEL = atoi (argv[1]);

  size(20.5);
  origin(0., -L0/2., -L0/2.);
  init_grid(64);

  mu = muv;

  cs.refine = cs.prolongation = fraction_refine; //cs is a volume fraction field

  TOLERANCE = 1e-4;
  DT=1e-3; //Time step
  
  run();
}

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*0.2/Re;
}

//BCs

//Inflow
u.n[left] = dirichlet(Uin*cs[]);
p[left] = neumann(0.);
pf[left] = neumann(0.);

//Outflow
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

event init (t = 0) {
  if (!restore (file = "restart")) {
    refine (fabs(y) < 0.6*Lyy && fabs(z) < 0.6*Lzz && level<LEVEL);

    //Channel
    vertex scalar phi[];
    foreach_vertex() {
      phi[] = intersection(Lyy/2.-y, Lyy/2.+y); //intersection(a,b) = min(a,b)
      phi[] = intersection(phi[], Lzz/2.-z);
      phi[] = intersection(phi[], Lzz/2.+z);
    }
    boundary ({phi});
    fractions (phi, cs);
    
    foreach() {
      u.x[] = cs[]; //Init velocity
      un[] = u.x[];
    }
  }
}

event velocity (i++) {
  foreach()
    foreach_dimension()
      u.x[] = cs[] * u.x[];
  boundary ((scalar *){u});
}

event logfile (i += 10; t <= 10) {
  double du = change (u.x, un);
  fprintf (ferr, "%d %g %g\n", i, t, du);
  fflush (ferr);
  if (i > 0 && du < 1e-6){ //Convergence criteria
    event("vtk_file");
    return 1; //Stop the simulation
  }
}

event profiles (t = end) {
  scalar * my_list = {u.x}; //list of scalars I want to export
  int len_my_list = list_len(my_list); //lenght of my list
  int np = 100;
  double v[(np+1)*len_my_list]; //number of interpolated points (np+1) times number of scalars (len_my_list)
  
  //line 1
  coord a[np+1];
  for (int n = 0; n <= np; n++) {
    a[n].x = 20.;
    a[n].y = 0.;
    a[n].z = -0.5 + (1./np) * n;
  }
  interpolate_array(my_list, a, np+1, v, true);

  if (pid()==0) {
    FILE * fp1 = fopen ("profile_1", "w");
    for (int n = 0; n <= np; n++) {
      fprintf(fp1, "%g %g %g %g\n", a[n].x, a[n].y, a[n].z,
	      v[n*len_my_list]);
    }
    fclose(fp1);
  }
}

//event movie (t += 0.1) {
//  view (fov = 22.4578, quat = {-0.707107,-0,-0,0.707107}, tx = -0.5, ty = 0.,
//  	bg = {0.3,0.4,0.6}, width = 600, height = 600, samples = 1);
//
//  clear();
//  squares("u.x", min=0, max=1.5, alpha = 0, n = {0,1,0});
//  save("movie_ux.mp4");
//}

event snapshot (t ++) {
  char name[80];
  sprintf(name, "dump-%06g",t);
  dump (file = name);
}

event vtk_file (t += 1){
  char subname[80]; sprintf(subname, "duct");

  output_vtu_MPI( (scalar *) {p, cs},
  (vector *) {u},
  (vector *) {uf}, subname, t + dt );
}