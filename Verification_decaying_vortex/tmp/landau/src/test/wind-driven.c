/**
# Wind-driven lake

This is a simple test case of a wind-driven lake where we can compare
results with an analytical solution. For the bottom of the domain we
impose *no-slip* condition (that is the default), for the top we
impose a Neumann condition (see [viscous friction between
layers](/src/multilayer.h#viscous-friction-between-layers)
for details). 

We run the test case for three different solvers: multilayer
Saint-Venant, layered hydrostatic, layered non-hydrostatic. */

#include "grid/multigrid1D.h"
#if !LAYERS
#  include "saint-venant.h"
#else
#  include "layered/hydro.h"
#  if NH
#    include "layered/nh.h"
#  endif
#  include "layered/remap.h"
#endif

/**
There are five parameters $L$, $h$, $g$, $\dot{u}$, $\nu$. Only three are
independent e.g.
$$
  a = \frac{L}{h},
$$
$$
  \text{Re} = \frac{\dot{u} h^2}{\nu} = s \frac{gh^3}{\nu^2},
$$
$$
  s = \frac{\dot{u} \nu}{gh} = \frac{\dot{u}^2 h}{\text{Re} g}
$$
where $a$ is the aspect ratio, Re the Reynolds number and $s$ the slope of the
free-surface. A characteristic time scale is
$$
t_{\nu} = \frac{h^2}{\nu}
$$
We choose a small Reynolds number and a small slope. */

double Re = 10.;
double s = 1./1000.;

double du0;
scalar duv[];

int main()
{
  L0 = 10.;
  X0 = -L0/2.;
  G = 9.81;
  N = 64;
  nu = sqrt(s*G/Re);
  du0 = sqrt(s*Re*G);
  dut = duv;

  /**
  We vary the number of layers. */
  
  for (nl = 4; nl <= 32; nl *= 2)
    run();
}

/**
We set the initial water level to 1 and set the surface stress. */

event init (i = 0) {
  foreach() {
#if !LAYERS
    h[] = 1.;
#else
    foreach_layer()
      h[] = 1./nl;
#endif
    duv[] = du0*(1. - pow(2.*x/L0,10));
  }
}

/**
We compute the error between the numerical solution and the analytical
solution. */

#define uan(z)  (du0*(z)/4.*(3.*(z) - 2.))

event error (t = 10./nu)
{
  int i = 0;
  foreach() {
    if (i++ == N/2) {
      double z = zb[], emax = 0.;
#if !LAYERS
      int l = 0;
      for (vector u in ul) {
	double e = fabs(u.x[] - uan (z + h[]*layer[l]/2.));
	if (e > emax) 
	  emax = e;
	z += h[]*layer[l++];
      }
#else
      foreach_layer() {
	double e = fabs(u.x[] - uan (z + h[]/2.));
	if (e > emax)
	  emax = e;
	z += h[];
      }
#endif
      fprintf (stderr, "%d %g\n", nl, emax);
    }
  }
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
#if !LAYERS
  fprintf (fp, "\n");
  foreach()
    fprintf (fp, "%g %g %g\n", x, zb[] + h[], zb[]);
#else
  int i = 4;
  for (scalar h in hl)
    fprintf (fp, ", '' u 1:%d w l lw 2 t ''", i++);
  fprintf (fp, "\n");
  foreach() {
    double H = 0.;
    foreach_layer()
      H += h[];
    fprintf (fp, "%g %g %g", x, zb[] + H, zb[]);
    double z = zb[];
    foreach_layer()
      fprintf (fp, " %g", z += h[]);
    fprintf (fp, "\n");
  }
#endif  
  fprintf (fp, "e\n\n");
  //  fprintf (fp, "pause 0.05\n");
  fflush (fp);
}
#endif

/**
We save the horizontal velocity profile at the center of the domain
and the two components of the velocity field for the case with 32
layers.*/

event output (t = end) {
  char name[80];
  sprintf (name, "uprof-%d", nl);
  FILE * fp = fopen (name, "w");
  int i = 0;
#if !NH && LAYERS
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
#endif // !NH && LAYERS
  foreach() {
    if (i++ == N/2) {
#if !LAYERS
      int l = 0;
      double z = zb[] + h[]*layer[l]/2.;
      for (vector u in ul)
	fprintf (fp, "%g %g\n", z, u.x[]), z += h[]*layer[l++];
#else
      double z = zb[];
      foreach_layer()
	fprintf (fp, "%g %g\n", z + h[]/2., u.x[]), z += h[];
#endif
    }
    if (nl == 32) {
      double z = zb[];
#if !LAYERS
      int l = 0;
      scalar w;
      vector u;
      for (w,u in wl,ul) {
	printf ("%g %g %g %g\n", x, z + h[]*layer[l]/2., u.x[], w[]);
	z += layer[l++]*h[];
      }
#else
      foreach_layer()
	printf ("%g %g %g %g\n", x, z + h[]/2., u.x[], w[]), z += h[];
#endif
      printf ("\n");
    }
  }
  fclose (fp);
#if !NH && LAYERS
  delete ({w});
#endif
}

/**
## Results

~~~gnuplot Numerical and analytical velocity profiles at the center of the lake.
set xr [0:1]
set xl 'z'
set yl 'u'
set key left top
G = 9.81
s = 1./1000.
Re = 10.
du0 = sqrt(s*Re*G)
plot [0:1]du0*x/4.*(3.*x-2.) t 'analytical', \
          'uprof-4' pt 5 t '4 layers', \
          'uprof-8' pt 6 t '8 layers', \
          'uprof-16' pt 9 t '16 layers', \
          'uprof-32' pt 10 t '32 layers'
~~~

~~~gnuplot Convergence of the error between the numerical and analytical solution with the number of layers.
reset
set cbrange [1:2]
set logscale
set xlabel 'Number of layers'
set ylabel 'max|e|'
set xtics 4,2,32
set grid
fit a*x+b 'log' u (log($1)):(log($2)) via a,b
plot [3:36]'log' u 1:2 pt 7 t '', \
     exp(b)*x**a t sprintf("%.2f/N^{%4.2f}", exp(b), -a) lt 1
~~~

~~~gnuplot Velocity field (32 layers).
reset
unset key
set xlabel 'x'
set ylabel 'z'
scale = 10.
plot [-5:5][0:1]'out' u 1:2:($3*scale):($4*scale) w vectors
~~~

## See also

* [Similar test with Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/lake.html#river)
*/
