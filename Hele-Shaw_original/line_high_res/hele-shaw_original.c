#include "hele-shaw.h"

#define MAXLEVEL 12

double mu1 = 2.5, mu2 = 1e-3, k = 0.1;
int iteration=0;
scalar f[];
scalar * tracers = {f};

p[left]  = dirichlet(1e-3);
p[right] = dirichlet(0);
f[left]  = 1.;

int main()
{
  size (7.5e-2);
  init_grid (1 << MAXLEVEL);
  run();
}

event init (i = 0) {
  foreach()
    f[] = (x < 1e-3)*(1. - 1e-3*(1. + noise()));
}

event logfile (i++)
{
  stats s = statsf (f);
  fprintf (stderr, "%d %g %d %g %g %g\n", 
	   i, t, mgp.i, s.sum, s.min, s.max);
}

event output (t += 10; t <= 50)
{
}

event coefficients (i++)
{
  foreach_face() {
    double ff = (f[] + f[-1])/2.;
    beta.x[] = - k/(mu1 + clamp(ff,0,1)*(mu2 - mu1));
  }
  boundary ((scalar *){beta});
}

#if TREE
event adapt (i++) {
  double tolerance = 0.02;
  adapt_wavelet ({f}, &tolerance, MAXLEVEL);
  event ("coefficients");
}
#endif


//Output
#include "output_fields/output_vtu_foreach.h"
event vtk_file (t += 0.1; t<=100)
{
  int nf = iteration;
  scalar l[];
  foreach()
  l[] = level;

  char name[80], subname[80];
  FILE *fp;
  sprintf(name, "hs_%4.4d_n%3.3d.vtu", nf, pid());
  fp = fopen(name, "w");

  output_vtu_bin_foreach((scalar *) {l, f, p}, (vector *) {u}, 64, fp, false);
  fclose(fp);
  @if _MPI
  if (pid() == 0) {
    sprintf(name, "hs_%4.4d.pvtu", nf);
    sprintf(subname, "hs_%4.4d", nf);
    fp = fopen(name, "w");
    output_pvtu_bin((scalar *) {l, f, p}, (vector *) {u}, 64, fp, subname);
    fclose(fp);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  @endif

        iteration++;
}

event snapshot (i += 1000) {
  char name[80];
  sprintf (name, "dump-%d", i);
  scalar pid[];
  foreach()
  pid[] = fmod(pid()*(npe() + 37), npe());
  boundary ({pid});
  dump (name);
}