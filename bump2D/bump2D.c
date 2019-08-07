//Bouncing Saint-Venant bump
#include "conservation.h"
scalar h[];
vector q[];
scalar * scalars = {h};
vector * vectors = {q};
#define LEVEL 7

double G = 1.;

void flux (const double * s, double * f, double e[2])
{
  double h = s[0], qx = s[1], u = qx/h, qy = s[2];
  f[0] = qx;
  f[1] = qx*u + G*h*h/2.;
  f[2] = qy*u;
  // min/max eigenvalues
  double c = sqrt(G*h);
  e[0] = u - c; // min
  e[1] = u + c; // max
}
event init (i = 0)
{
  foreach()
    h[] = 0.1 + 1.*exp(-200.*(x*x + y*y));
}

event logfile (i++) {
  stats s = statsf (h);
  fprintf (stderr, "%g %d %g %g %.8f\n", t, i, s.min, s.max, s.sum);
}

event graphs (i++) {
  stats s = statsf (h);
  fprintf (stderr, "%g %g %g\n", t, s.min, s.max);
}

event end (t = 4) {
  printf ("i = %d t = %g\n", i, t);
}

event outputfile (t+= 4./300.) {

    output_ppm (h, linear = true);
    scalar l[];
    foreach()
    l[] = level;
    static FILE * fp = fopen ("grid.ppm", "w");
    output_ppm (l, fp, min = 0, max = LEVEL);
  /* check symmetry */
  foreach() {
    double h0 = h[];
    point = locate (-x, -y);
    //    printf ("%g %g %g %g %g\n", x, y, h0, h[], h0 - h[]);
    assert (fabs(h0 - h[]) < 1e-12);
    point = locate (-x, y);
    assert (fabs(h0 - h[]) < 1e-12);
    point = locate (x, -y);
    assert (fabs(h0 - h[]) < 1e-12);
  }
}

#if TREE
event adapt (i++) {
  astats s = adapt_wavelet ({h}, (double[]){1e-3}, LEVEL);
  fprintf (stderr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
}
#endif



int main()
{
    origin (-0.5, -0.5);
    init_grid (1 << LEVEL);
    run();
}