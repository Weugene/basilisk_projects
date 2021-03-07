#include "hele-shaw.h"
#include "embed.h"
#define MAXLEVEL 8
#define BOOL_CYLINDER (sq(x-L0/2.) + sq(y-L0/2.) < sq(1e-2) )
#define CYLINDER (0.25e-2 - sqrt(sq(x - L0/2.) + sq(y - L0/2.)))

double mu1 = 2.5, mu2 = 1e-3, k = 0.1;
int iteration=0;
scalar f[];
scalar * tracers = {f};
face vector muc[];


int main()
{
  size (7.5e-2);
  init_grid (1 << MAXLEVEL);
  run();
}

p[right] = dirichlet(0);
p[left] = dirichlet(0);
p[top] = dirichlet(0);
p[bottom] = dirichlet(0);

p[embed] = dirichlet(1e-3);

vertex scalar eps[];

event init (i = 0) {
    int it = 0;
    do {
        it++;
        foreach()
        {
            f[] = BOOL_CYLINDER * (1. - 1e-3 * (1. + noise()));
            p[] = 1.e-3 * f[];
        }
        foreach_vertex()
        {
            eps[] = -CYLINDER;
        }
        boundary({eps});
        fractions(eps, cs, fs); //eps<0 in cyl, eps>0 out of cyl
        if (it>=10) printf("WARNING: does not converge... ");
    }while (adapt_wavelet({cs, f}, (double []){0.001, 0.001},
                          maxlevel = MAXLEVEL, minlevel = 3).nf != 0 && it <= 10);

    p[embed] = dirichlet (1e-3);
    boundary (all);
    event("coefficients");
    event ("vtk_file");
}

event logfile (i++)
{
  stats s = statsf (f);
  fprintf (stderr, "%d %g %d %g %g %g\n", 
	   i, t, mgp.i, s.sum, s.min, s.max);
}

event coefficients (i++)
{
    double ff;
  foreach_face() {
    ff = (f[] + f[-1])/2.;
    muc.x[] = (mu1 + clamp(ff,0,1)*(mu2 - mu1));
    beta.x[] = -k / muc.x[]*fs.x[];
  }
  boundary ((scalar *){beta});
}
#if TREE
event adapt (i++; t <= 50) {
  double tolerance = 0.02;
  adapt_wavelet ({f, cs}, (double[]){0.001, 0.001}, MAXLEVEL);
//  event ("coefficients");
}
#endif


//Output
#include "output_fields/output_vtu_foreach.h"
event vtk_file (t += 0.01; t<=100)
{
  int nf = iteration;
  scalar l[];
  foreach()
  l[] = level;

  char name[80], subname[80];
  FILE *fp;
  sprintf(name, "hs_%4.4d_n%3.3d.vtu", nf, pid());
  fp = fopen(name, "w");

  output_vtu_bin_foreach((scalar *) {l, f, p, cs, eps}, (vector *) {u, fs, muc, beta}, 64, fp, false);
  fclose(fp);
  @if _MPI
  if (pid() == 0) {
    sprintf(name, "hs_%4.4d.pvtu", nf);
    sprintf(subname, "hs_%4.4d", nf);
    fp = fopen(name, "w");
    output_pvtu_bin((scalar *) {l, f, p, cs, eps}, (vector *) {u, fs, muc, beta}, 64, fp, subname);
    fclose(fp);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  @endif

        iteration++;
}
