#include "grid/cartesian1D.h"
#include "conservation.h"

scalar h[];
vector q[];
scalar * scalars = {h};
vector * vectors = {q};

double G = 1.;

void flux (const double * s, double * f, double e[2])
{
  double h = s[0], q = s[1], u = q/h;
  f[0] = q;
  f[1] = q*u + G*h*h/2.;
  f[2] = s[2]*u;
  // min/max eigenvalues
  double c = sqrt(G*h);
  e[0] = u - c; // min
  e[1] = u + c; // max
}

int main()
{
  theta = 1.;
  origin (-0.5, -0.5);
  init_grid (500);
  run();
}

event init (i = 0)
{
  foreach()
    h[] = 0.1 + exp(-200.*x*x);
}

event logfile (t += 0.1; t <= 0.7) {
  foreach()
    fprintf (stderr, "%g %g %.6f\n", x, h[], q.x[]);
  fprintf (stderr, "\n");
}
