/* This is similar to gerris/test/poisson/circle */

#include "utils.h"
#include "poisson.h"

scalar a[], b[], res[], dp[];

double solution (double x, double y, double z)
{
  return cos(3.*pi*x)*cos(3.*pi*y)*cos(3.*pi*z);
}

void solve (int depth)
{
  /* Dirichlet condition on all boundaries */
  foreach_dimension() {
    a[right] = dirichlet (solution(x, y, z));
    a[left]  = dirichlet (solution(x, y, z));
  }
  /* homogeneous conditions for dp */
  foreach_dimension() {
    dp[right] = dirichlet(0);
    dp[left]  = dirichlet(0);
  }
  
  origin (-0.5, -0.5, -0.5);
  int nrelax = 4;
  init_grid(1);

  refine (level < depth - 2 || level <= depth*(1. - sqrt(x*x + y*y + z*z)));
  
  foreach() {
    b[] = - 9.*dimension*pi*pi*cos(3.*pi*x)*cos(3.*pi*y)*cos(3.*pi*z);
    a[] = 0.;
  }
  boundary ({a});

  #define NITER 15
  clock_t start = clock(), iter[NITER];
  double maxres[NITER];
  const scalar lambda[] = 0.;
  struct Poisson p;
  p.a = a; p.b = b; p.alpha = unityf; p.lambda = lambda;
  scalar * lres = {res};
  residual ({a}, {b}, lres, &p);
  for (int i = 0; i < NITER; i++) {
    mg_cycle ({a}, lres, {dp}, relax, &p, nrelax, 0, depth());
    residual ({a}, {b}, lres, &p);
    double max = 0.;
    foreach(reduction(max:max))
      if (fabs(res[]) > max)
	max = fabs(res[]);
    iter[i] = clock();
    maxres[i] = max;
  }
  for (int i = 0; i < NITER; i++) {
    fprintf (stderr, "residual %d %d %.1g\n", depth, i, maxres[i]);
    printf ("speed %d %d %g %g\n", depth, i, 
	    (iter[i] - start)/(double)CLOCKS_PER_SEC, maxres[i]);
  }

  double max = 0;
  foreach(reduction(max:max)) {
    double e = a[] - solution(x, y, z);
    if (fabs(e) > max) max = fabs(e);
    //    printf ("%g %g %g %g %g %g\n", x, y, a[], b[], res[], e);
  }
  fprintf (stderr, "max error %d %g\n", depth, max);
}

int main (int argc, char ** argv)
{
#if dimension <= 2
  for (int depth = 7; depth <= 10; depth++)
#else
  for (int depth = 5; depth <= 7; depth++)
#endif
    solve (depth);
}
