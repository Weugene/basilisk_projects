// Bumps uniformly covered with water create ponds and streams (i.e. a
// lot of wetting/drying going on)

// use KDT database rather than analytical function
#define KDT 1

#if IMPLICIT
# include "saint-venant-implicit.h"
#else
# include "saint-venant.h"
#endif

#if KDT
# include "terrain.h"
#endif

#define LEVEL 7

int main()
{
  init_grid (1 << LEVEL);
  size (1000.);
  G = 9.81;
#if IMPLICIT
  DT = 10.;
#endif
  run();
}

#define zb(x,y) ((cos(pi*x/L0)*cos(pi*y/L0) + \
		  cos(3.*pi*x/L0)*cos(3.*pi*y/L0)) - 2.*x/1000.)

#if !KDT
void refine_zb (Point point, scalar zb)
{
  foreach_child()
    zb[] = zb(x,y);
}
#endif

event init (i = 0)
{
#if !KDT
  zb.refine = refine_zb; // updates terrain
  foreach() {
    zb[] = zb(x,y);
    h[] = 0.1;
  }
#else
  FILE * fp = popen ("xyz2kdt ponds", "w");
  for (double x = 0.; x <= 1000; x += 1.)
    for (double y = 0.; y <= 1000.; y += 1.)
      fprintf (fp, "%g %g %g\n", x, y, zb(x,y));
  pclose (fp);
  terrain (zb, "ponds", NULL);
  foreach()
    h[] = 0.1;
#endif
}

event friction (i++) {
#if IMPLICIT
  // quadratic bottom friction, coefficient 1e-4 (dimensionless)
  foreach() {
    double a = h[] < dry ? HUGE : 1. + 1e-4*dt*norm(q)/sq(h[]);
    foreach_dimension()
      q.x[] /= a;
  }
  boundary ((scalar *){q});
#else
  // quadratic bottom friction, coefficient 1e-4 (dimensionless)
  foreach() {
    double a = h[] < dry ? HUGE : 1. + 1e-4*dt*norm(u)/h[];
    foreach_dimension()
      u.x[] /= a;
  }
  boundary ((scalar *){u});
#endif
}

event logfile (i += 10) {
  stats s = statsf (h);
#if IMPLICIT
  scalar u[];
  foreach()
    u[] = h[] > dry ? q.x[]/h[] : 0.;
  norm n = normf (u);
#else
  norm n = normf (u.x);
#endif
  if (i == 0)
    fprintf (stderr, "t i h.min h.max h.sum u.x.rms u.x.max dt\n");
  fprintf (stderr, "%g %d %.4g %g %.8f %g %g %g %g\n", 
	   t, i, s.min, s.max, s.sum, dt, n.rms, n.max, dt);
#if !IMPLICIT
  assert (s.min > 0.);
#endif
}

event outputfile (t <= 1200.; t += 1200./8) {
#if !IMPLICIT  
  static int nf = 0;
  printf ("file: eta-%d\n", nf);
  output_field ({h, zb, u}, stdout, N, linear = true);

  scalar l[];
  foreach()
    l[] = level;
  printf ("file: level-%d\n", nf);
  output_field ({h, zb, l}, stdout, N);

  nf++;
#endif
}

// int event (t += 100)
//  output_matrix (h, stdout, N, true);

event adapt (i++) {
  // we dot this so that wavelets use the default bilinear
  // interpolation this is less noisy than the linear + gradient
  // limiters used in Saint-Venant not sure whether this is better
  // though.
  scalar h1[];
  foreach()
    h1[] = h[];
  boundary ({h1});
  adapt_wavelet ({h1}, (double[]){1e-2}, LEVEL, 4); 
}
