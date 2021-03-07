/**
# Sessile drop

A sessile drop is a drop of liquid at rest on a solid surface. In the
absence of gravity, the shape of the drop is controlled by surface
tension only. An important parameter is the "contact angle" $\theta$ between
the solid surface and the interface. In the absence of gravity, the
drop is hemispherical and it is easy to show that the relation between
the radius of the drop $R$ and its volume $V$ is (for two-dimensional
drops)
$$
V = R^2 (\theta - \sin\theta\cos\theta)
$$

To test this relation, a drop is initialised as a half-disk (i.e. the
initial contact angle is 90$^\circ$) and the contact angle is varied
between 15$^\circ$ and 165$^\circ$. The drop oscillates and eventually relaxes
to its equilibrium position. This equilibrium is exact to within
machine accuracy. The curvature along the interface is constant.

Note that shallower angles are [not accessible yet](/src/contact.h).

~~~gnuplot Equilibrium shapes for $15^\circ \leq \theta \leq 165^\circ$
set term push
set term @SVG size 640,180
set size ratio -1
unset key
unset xtics
unset ytics
unset border
plot 'out' w l, '' u (-$1):2 w l lt 1, 0 lt -1
set term pop
~~~
*/

//#include "grid/multigrid.h"
//#include "grid/quadtree.h"
//#include "navier-stokes/centered.h"
#include "../src_local/centered-weugene.h"
#include "contact.h"
#include "vof.h"
#include "tension.h"
#define MAXLEVEL 8
scalar f[], * interfaces = {f};

/**
To set the contact angle, we allocate a [height-function
field](/src/heights.h) and set the contact angle boundary condition on
its tangential component. */

vector h[];
double theta0 = 30;
h.t[bottom] = contact_angle (theta0*pi/180.);
scalar kappa_cur[];
int main()
{
  size (2);

  /**
  We use a constant viscosity. */

  const face vector muc[] = {.1,.1};
  mu = muc;

  /**
  We must associate the height function field with the VOF tracer, so
  that it is used by the relevant functions (curvature calculation in
  particular). */

  f.height = h;

  /**
  We set the surface tension coefficient and run for the range of
  contact angles. */

  f.sigma = 1.;

  for (theta0 = 15; theta0 <= 165; theta0 += 15)
    run();
}

/**
The initial drop is a quarter of a circle. */

event init (t = 0)
{
  fraction (f, - (sq(x) + sq(y) - sq(0.5)));
}

#if 0
event logfile (i++)
{
  fprintf (fout, "%g %g\n", t, normf(u.x).max);
}

event snapshots (t += 1)
{
  p.nodump = false;
  dump();
}
#endif

//Output
#include "../src_local/output_vtu_foreach.h"
event vtk_file (i += 100)
{
    char subname[80]; sprintf(subname, "ses");
    curvature (f, kappa_cur);
    output_vtu_MPI( (scalar *) {f, p, pf, kappa_cur}, (vector *) {a, u, uf}, subname);
}

/**
At equilibrium (t = 10 seems sufficient), we output the interface
shape and compute the (constant) curvature. */

event end (t = 10)
{
  output_facets (f, stdout);


  curvature (f, kappa_cur);
  stats s = statsf (kappa_cur);
  double R = s.volume/s.sum, V = 2.*statsf(f).sum;
  fprintf (ferr, "%d %g %.5g %.3g\n", N, theta0, R/sqrt(V/pi), s.stddev);
}



#if TREE
event adapt (i++) {
#if 1
  scalar f1[];
  foreach()
    f1[] = f[];
  boundary ({f1});
  adapt_wavelet ({f1}, (double[]){1e-3}, minlevel = 3, maxlevel = MAXLEVEL);
#else
  adapt_wavelet ({f}, (double[]){1e-4}, minlevel = 3, maxlevel = MAXLEVEL);
#endif
}
#endif

/**
We compare $R/R_0$ to the analytical expression, with $R_0=\sqrt{V/\pi}$.

~~~gnuplot
reset
set xlabel 'Contact angle (degrees)'
set ylabel 'R/R_0'
set arrow from 15,1 to 165,1 nohead dt 2
set xtics 15,15,165
plot 1./sqrt(x/180. - sin(x*pi/180.)*cos(x*pi/180.)/pi) t 'analytical', \
  'log' u 2:3 pt 7 t 'numerical'
~~~

## See also

* [Similar test with
   Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/sessile.html)
*/