/**
# Moving embed test case. */

# include "grid/quadtree.h"
# include "embed.h"
#include "output_vtu_foreach.h"
/**
## Mesh parameters */

# define LEVEL 9 
# define MAXLEVEL 9

/**
## Physical parameters */

# define Re 100 // Reynolds number
# define diamsolid 0.1 // Particle diameter
# define xinit (7.*diamsolid)
# define yinit (L0/2.)
# define vc 1. // Velocity of the cylinder
# define rhoval 1. // Fluid density
# define tval (vc*diamsolid/Re) // Viscosity

# define NDT 1.
# define T_CONV (diamsolid/vc)
# define T_DIFF (Re*(T_CONV))
# define T_MESH ((L0/(1 << MAXLEVEL))/(vc))

# include "navier-stokes/centered.h"
# include "view.h"

double xpos, ypos;

/**
## Particle shape */

void particle (scalar cs, face vector fs, double xp, double yp) {

  vertex scalar phi[];

  foreach_vertex()
    phi[] = sq (x - xp) + sq (y - yp) - sq (diamsolid/2.);

  boundary ({phi});
  fractions (phi, cs, fs);
  boundary ((scalar *) {cs, fs});
}

/**
## Code */

int main() {

  L0 = 1.;
  size (L0);
  origin (0., 0.);

  stokes = false;
  TOLERANCE = 1e-5;
  DT = min (T_MESH, min (T_CONV, T_DIFF))/NDT;
  
  N = 1 << MAXLEVEL;
  init_grid (N);
  run ();
}

/**
#### Boundary conditions */

u.n[left] = dirichlet (-vc);
u.t[left] = dirichlet (0);

u.n[right] = dirichlet (-vc);
u.t[right] = dirichlet (0);

/**
Embedded boundary condition. */

u.n[embed] = dirichlet (0.);
u.t[embed] = dirichlet (0.);
p[embed] = neumann (0.);

/**
#### Initial conditions */

face vector muv[];

event init (t = 0) {

  foreach() {
    p[] = 0.;
    u.x[] = -vc;
    u.y[] = 0.;
  }

  boundary ((scalar *) {u, p});

  /**
  We define the initial embedded geometry. */

  xpos = (xinit);
  ypos = (yinit);
  particle (cs, fs, xpos, ypos);
  
  fractions_cleanup (cs, fs);
  foreach_face()
    if (uf.x[] && !fs.x[])
      uf.x[] = 0.;
  boundary ((scalar *){cs, fs, uf});
    
  /**
  Viscosity. */
  
  mu = muv;
}

/**
#### Properties */

event properties (i++) {

  foreach_face()
    muv.x[] = fm.x[]*tval;
  boundary ((scalar *) {muv});
}

/**
## Visualization */

event drag_coeff (i++) {

  coord Fp, Fmu;
  embed_force (p, u, mu, &Fp, &Fmu);

  double Sproj = pi/4.*sq (diamsolid);
  double Ep = 0.5*rhoval*sq ((vc))*Sproj;
  double Cd = (Fp.x + Fmu.x)/Ep;
  double Cl = (Fp.y + Fmu.y)/Ep;

  static FILE * fp = fopen ("coefficients_eb.dat", "w");
  fprintf (fp, "%d %g %g %g\n", i, t, Cd, Cl);
  fprintf (ferr, "coeffs: %d %g %g %g\n", i, t, Cd, Cl);
  fflush (fp);
}

//event movie (i++) {
//
//  scalar omega[];
//  vorticity (u, omega);
//  fprintf(ferr, "1");
//  view (fov = 3.6, quat = {0,0,0,1},
//	tx = -xpos/L0, ty = -0.5, bg = {1,1,1},
//	width = 600, height = 600, samples = 1);
//fprintf(ferr, "2");
//  cells ();
//fprintf(ferr, "3");
//  draw_vof ("cs", "fs", lw = 5);
//  squares ("omega", min = -10, max = 10, linear = false);
//  save ("cells-omega.mp4");
//
//  view (fov = 3.6, quat = {0,0,0,1},
//	tx = -xpos/L0, ty = -0.5, bg = {1,1,1},
//	width = 600, height = 600, samples = 1);
//  cells ();
//  draw_vof ("cs", "fs", lw = 5);
//  squares ("u.x", linear = false);
//  save ("cells-ux.mp4");
//
//  view (fov = 3.6, quat = {0,0,0,1},
//	tx = -xpos/L0, ty = -0.5, bg = {1,1,1},
//	width = 600, height = 600, samples = 1);
//  cells ();
//  draw_vof ("cs", "fs", lw = 5);
//  squares ("u.y", linear = false);
//  save ("cells-uy.mp4");
//
//  view (fov = 3.6, quat = {0,0,0,1},
//	tx = -xpos/L0, ty = -0.5, bg = {1,1,1},
//	width = 600, height = 600, samples = 1);
//  cells ();
//  draw_vof ("cs", "fs", lw = 5);
//  squares ("p", linear = false);
//  save ("cells-p.mp4");
//
//  view (fov = 23., quat = {0,0,0,1},
//	tx = -0.5, ty = -0.5, bg = {1,1,1},
//	width = 600, height = 600, samples = 1);
//  draw_vof ("cs", "fs", lw = 5);
//  squares ("omega", linear = false);
//  save ("omega.mp4");
//
//  view (fov = 23., quat = {0,0,0,1},
//	tx = -0.5, ty = -0.5, bg = {1,1,1},
//	width = 600, height = 600, samples = 1);
//  draw_vof ("cs", "fs", lw = 5);
//  squares ("u.x", linear = false);
//  save ("ux.mp4");
//
//  view (fov = 23., quat = {0,0,0,1},
//	tx = -0.5, ty = -0.5, bg = {1,1,1},
//	width = 600, height = 600, samples = 1);
//  draw_vof ("cs", "fs", lw = 5);
//  squares ("u.y", linear = false);
//  save ("uy.mp4");
//
//  view (fov = 23., quat = {0,0,0,1},
//	tx = -0.5, ty = -0.5, bg = {1,1,1},
//	width = 600, height = 600, samples = 1);
//  draw_vof ("cs", "fs", lw = 5);
//  squares ("p", linear = false);
//  save ("p.mp4");
//}


//event vtk_file (i+=10)
event vtk_file (t += 0.1)
{
char subname[80]; sprintf(subname, "fixed_eb");
scalar l[]; foreach() l[] = level;

fprintf(ferr, "output_vtu_MPI\n");
output_vtu_MPI(subname, (iter_fp) ? t + dt : 0, list = (scalar *) {p, fs, l}, vlist = (vector *) {u});
}

event stop_run (t = 20) {
  return 1.;
}

/**
## Plots

![Vorticity *omega*](test-embed-fixed-cylinder-2D-Re-100/omega.mp4)
![Velocity *u.x*](test-embed-fixed-cylinder-2D-Re-100/ux.mp4)
![Velocity *u.y*](test-embed-fixed-cylinder-2D-Re-100/uy.mp4)
![Pressure *p*](test-embed-fixed-cylinder-2D-Re-100/p.mp4)
*/
