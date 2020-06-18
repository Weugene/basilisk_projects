/**
# Transcritical flow over a bump with multiple layers

This test case is similar to [transcritical.c]() but using simpler
(and more physical) boundary conditions. Both hydrostatic and
non-hydrostatic solutions are obtained. The non-hydrostatic solution
is compared with [a DNS solution](/src/examples/gaussian-ns.c) using
the Navier-Stokes/VOF solver.

More details are given in [Popinet (2019)](/Bibliography#popinet2019). */

#include "grid/multigrid1D.h"
#if !LAYERS
# include "saint-venant.h"
#else // LAYERS
# include "layered/hydro.h"
# if !HYDRO
#   include "layered/nh.h"
# endif
# include "layered/remap.h"
# include "layered/perfs.h"
#endif // LAYERS

/**
The primary parameters are the flow rate, water level, viscosity and bump
amplitude. */

#define QL 1.
#define HR 0.6
#define NU 1e-2
#define BA 0.4

int main()
{

  /**
  The domain is 30 meters long. */
  
  L0 = 30.;
  G = 9.81;
  N = 512; // less damping with 1024

  /**
  The viscosity is set to NU and 20 layers are used to discretise
  vertically. */
  
  nu = NU;
  nl = 20; // going to 30 changes very little
  run();
}

/**
## Initialisation and boundary conditions

We create a field *hc* to check convergence on *h*. */

scalar hc[];

/**
The inflow is a parabolic profile with a total flow rate *Q*. The
function below computes the height *zc* of the middle of the layer and
returns the corresponding velocity. */

#if !LAYERS
double uleft (Point point, scalar s, double Q)
{
  double zc = zb[];
  for (int l = 0; l < s.l; l++)
    zc += layer[l]*h[];
  zc += layer[s.l]*h[]/2.;
  return 3./2.*Q/h[]*(1. - sq(zc/h[] - 1.));
}
#else

/**
For the multilayer solver, we need some trickery to define the inflow
velocity profile as a function of $z$. We need to sum the layer
thicknesses to get the total depth `H` and the height `zc` of the
current layer (of index `point.l`). */

double uleft (Point point, scalar s, double Q)
{
  double H = 0.;
  double zc = zb[];
  for (int l = - point.l; l < nl - point.l; l++) {
    H += h[0,0,l];
    if (l < 0)
      zc += h[0,0,l];
  }
  zc += h[]/2.;
  return 3./2.*Q/H*(1. - sq(zc/H - 1.));
}
#endif

/**
We initialise the topography and the initial thickness of each layer *h*. */

event init (i = 0) {
  foreach() {
    zb[] = BA*exp(- sq(x - 10.)/5.);
    hc[] = HR - zb[];
#if !LAYERS
    h[] = hc[];
#else
    foreach_layer()
      h[] = hc[]/nl;
#endif
  }

  /**
  The height of the free-surface $\eta$ is imposed on the right
  boundary. */
  
  eta[right] = dirichlet(HR);

  /**
  Boundary conditions are set for the inflow velocity and the outflow
  of each layer. */
  
#if !LAYERS
  for (vector u in ul) {
    u.n[left] = dirichlet (uleft (point, _s, QL*(t < 10. ? t/10. : 1.)));
    u.n[right] = neumann(0.);
  }
  h[right] = dirichlet(HR);
#else
  u.n[left] = dirichlet (uleft (point, _s, QL*(t < 10. ? t/10. : 1.)));
  u.n[right] = neumann(0.);
  h[right] = HR/nl;
#endif

  /**
  In the non-hydrostatic case, a boundary condition is required for
  the non-hydrostatic pressure $\phi$ of each layer. */
  
#if !HYDRO && LAYERS
  phi[right] = dirichlet(0.);
#endif
}

/**
We can optionally add horizontal viscosity. */

#if 0
event viscous_term (i++)
{
  // add horizontal viscosity (small influence)
  foreach_layer() {
    scalar d2u[];
    foreach()
      d2u[] = (u.x[1] + u.x[-1] - 2.*u.x[])/sq(Delta);
    foreach()
      u.x[] += dt*nu*d2u[];
    foreach()
      d2u[] = (w[1] + w[-1] - 2.*w[])/sq(Delta);
    foreach()
      w[] += dt*nu*d2u[];
  }
  boundary ({u, w});
}
#endif

/**
We check for convergence. */

event logfile (t += 0.1; i <= 100000) {
#if !LAYERS
  double dh = change (h, hc);
#else
  scalar H[];
  foreach() {
    H[] = 0.;
    foreach_layer()
      H[] += h[];
  }
  double dh = change (H, hc);
#endif
  if (i > 0 && dh < 1e-5)
    return 1;
}

/**
Uncomment this part if you want on-the-fly animation. */

#if 0
event gnuplot (i += 20) {
  static FILE * fp = popen ("gnuplot 2> /dev/null", "w");
  if (i == 0)
    fprintf (fp, "set term x11\n");
  fprintf (fp,
	   "set title 'nl = %d, t = %.2f'\n"
	   "p [%g:%g][0:]'-' u 1:3:2 w filledcu lc 3 t '',"
	   " '' u 1:(-1):3 t '' w filledcu lc -1", nl, t,
	   X0, X0 + L0);
  int i = 4;
  foreach_layer()
    fprintf (fp, ", '' u 1:%d w l lw 2 t ''", i++);
  fprintf (fp, "\n");
  foreach_leaf() {
    double H = 0.;
    foreach_layer()
      H += h[];
    fprintf (fp, "%g %g %g", x, zb[] + H, zb[]);
    double z = zb[];
    foreach_layer() {
      fprintf (fp, " %g", z);
      z += h[];
    }
    fprintf (fp, "\n");
  }
  fprintf (fp, "e\n\n");
  //  fprintf (fp, "pause 0.05\n");
  fflush (fp);
}
#endif

/**
## Outputs

At the end of the simulation we save the profiles. */

event profiles (t += 5)
{
  foreach_leaf() {
#if !LAYERS
    double H = h[];
#else
    double H = 0.;
    foreach_layer()
      H += h[];
#endif
    fprintf (stderr, "%g %g %g\n", x, zb[] + H, zb[]);
  }
  fprintf (stderr, "\n\n");
}

/**
We also save the velocity and non-hydrostatic pressure fields. */

event output (t = end)
{
#if HYDRO  
  scalar w = new scalar[nl];
  vertical_velocity (w);
  foreach() {
    double wm = 0.;
    foreach_layer() {
      double w1 = w[];
      w[] = (w1 + wm)/2.;
      wm = w1;
    }
  }
  boundary ({w});
#endif // HYDRO
  foreach_leaf() {
    double z = zb[];
#if HYDRO
    printf ("%g %g %g %g\n", x, z, u.x[], w[]);
    foreach_layer() {
      z += h[];
      printf ("%g %g %g %g\n", x, z, u.x[], w[]);
    }
#elif !LAYERS
    vector u = ul[0];
    scalar w = wl[0];
    printf ("%g %g %g %g\n", x, z, u.x[], w[]);
    int l = 0;
    for (u,w in ul,wl) {
      z += layer[l++]*h[];
      printf ("%g %g %g %g\n", x, z, u.x[], w[]);
    }
#else // LAYERS
    printf ("%g %g %g %g %g\n", x, z, u.x[], w[], phi[]);
    foreach_layer() {
      z += h[];
      printf ("%g %g %g %g %g\n", x, z, u.x[], w[], phi[]);
    }
#endif // LAYERS   
    printf ("\n");
  }
#if HYDRO
  delete ({w});
#endif
  printf ("# end = %g\n", t);
}

/**
## Results

We compare the hydrostatic results obtained with the [layered
solver](/src/layered/hydro.h) and those obtained with the
[Saint-Venant multilayer code](/src/multilayer.h).

~~~gnuplot Evolution of the free surface (hydrostatic).
set yr [0.5:1.2]
set xlabel 'x'
set ylabel 'z'
plot '../gaussian-hydro/log' u 1:2 w l t 'Multilayer', \
     '../gaussian-stvt/log' u 1:2 w l t 'De Vita et al.'
~~~

The results which follow are obtained with the non-hydrostatic solver.

~~~gnuplot Evolution of the free surface (non-hydrostatic)
plot 'log' u 1:2 w l t ''
~~~

We compare the results obtained with the layered solver to those
obtained with the Navier-Stokes VOF solver (and without metric).

~~~gnuplot Evolution of the free surface for both solvers
plot 'log' w l t 'Multilayer', \
     '../../examples/gaussian-ns/log' w l t 'Navier-Stokes VOF'
~~~

~~~gnuplot Final free surface profiles
plot 'log' index 13 w l t 'Multilayer', \
     '../../examples/gaussian-ns/log' index 'prof70' w l t 'Navier-Stokes VOF', \
     'gaussian.nometric' w l t 'no metric'
~~~

~~~gnuplot Horizontal velocity field { width=100% }
set term PNG enhanced font ",15" size 2048,512
set yr [0:1.2]
set output 'vel.png'
set pm3d
set pm3d map interpolate 10,4
unset key
# jet colormap
set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1,	\
0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625 1 0.9333 0, 0.75 1 0.4392 0,	\
0.875 0.9333 0 0, 1 0.498 0 0 )
set contour base
set cntrparam levels discrete 0
set cntrlabel onecolor
set cntrparam bspline
splot 'out' u 1:2:3 lt 3 lc rgb "#ffffff" lw 2
~~~

~~~gnuplot Vertical velocity field { width=100% }
set term PNG enhanced font ",15" size 2048,512
set output 'w.png'
# set cbrange [-0.1:0.1]
unset contour
splot 'out' u 1:2:4
~~~

~~~gnuplot Non-hydrostatic pressure { width=100% }
set term PNG enhanced font ",15" size 2048,512
set output 'phi.png'
# set cbrange [-0.2:0.35]
splot 'out' u 1:2:5
~~~

## See also

* [Same case with the Navier-Stokes/VOF solver](/src/examples/gaussian-ns.c)
*/
