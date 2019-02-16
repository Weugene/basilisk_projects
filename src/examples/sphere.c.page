/**
# Vortex shedding behind a sphere at Reynolds = 300

We solve the Navier--Stokes equations on an adaptive octree. */

#include "grid/octree.h"
#include "navier-stokes/centered.h"

/**
We will use the $\lambda_2$ criterion of [Jeong and Hussain,
1995](/src/references.bib#jeong1995) for vortex detection. */

#include "lambda2.h"

/**
This is the maximum level of refinement i.e. an equivalent maximum
resolution of $256^3$. */

int maxlevel = 8;

/**
The domain size is $16^3$. We move the origin so that the center of
the unit sphere is not too close to boundaries. */

int main() {
  L0 = 16.;
  origin (-2, -L0/2., -L0/2.);

  /**
  The viscosity is just $1/Re$, because we chose a sphere of diameter
  unity and an unit inflow velocity. */
  
  const face vector muc[] = {1./300,1./300,1./300};
  mu = muc;

  init_grid (64);
  run();
}

/**
The boundary conditions are inflow with unit velocity on the
left-hand-side and outflow on the right-hand-side. */

u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

/**
We define a new, no-slip boundary condition (i.e. Dirichlet condition
for the tangential velocity *u.t*) for the sphere. */

bid sphere;
u.t[sphere] = dirichlet(0.);
u.r[sphere] = dirichlet(0.);

event init (t = 0) {

  /**
  We initially refine only in a sphere, slightly larger than the solid
  sphere. */

  refine (x*x + y*y + z*z < sq(0.6) && level < maxlevel);

  /**
  We define the unit sphere by masking. */
  
  mask (x*x + y*y + z*z < sq(0.5) ? sphere : none);

  /**
  We set the initially horizontal velocity to unity everywhere. */
  
  foreach()
    u.x[] = 1.;
}

/**
We log the number of iterations of the multigrid solver for pressure
and viscosity. */

event logfile (i++)
  fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);

/**
We use *gfsview* and *ppm2gif* to create the animated isosurface of
$\lambda_2$ for $30 <= t <= 60$. */

event movies (t = 30; t += 0.25; t <= 60) {
  static FILE * fp =
    popen ("gfsview-batch3D sphere.gfv | ppm2gif > sphere.gif", "w");

  /**
  Here we compute two new fields, $\lambda_2$ and the vorticity
  component in the $y-z$ plane. */
  
  scalar l2[], vyz[];
  foreach()
    vyz[] = ((u.y[0,0,1] - u.y[0,0,-1]) - (u.z[0,1] - u.z[0,-1]))/(2.*Delta);
  boundary ({vyz});
  lambda2 (u, l2);
  output_gfs (fp);

  /**
  We tell *gfsview-batch3D* to save the images. */
  
  fprintf (fp, "Save stdout { format = PPM width = 640 height = 480}\n");
}

/**
![Animation of the $\lambda_2$ vortices coloured with the vorticity
 component aligned with the flow.](sphere/sphere.gif)

We set an adaptation criterion with an error threshold of 0.02 on all
velocity components. */

event adapt (i++) {
  astats s = adapt_wavelet ((scalar *){u}, (double[]){0.02,0.02,0.02},
			    maxlevel, 4);
  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
