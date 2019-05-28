#include "hele-shaw.h"
#define MAXLEVEL 10
#define RAD 1.e-2
#define RADINLET 0.5e-2
#define BOOL_CYLINDER (sq(x) + sq(y) < sq(RAD))
#define INLET      (sq(x) + sq(y) < sq(RADINLET))

double mu1 = 2.5, mu2 = 2.5, k = 1.;
int iteration=0;
scalar f[];
scalar * tracers = {f};

bid inletflow;
//p[inletflow] = dirichlet(0.);
//f[inletflow] = dirichlet(1.);

p[right] = dirichlet(0);
p[left] = dirichlet(0);
p[top] = dirichlet(0);
p[bottom] = dirichlet(0);

int main()
{
  size (7.5e-2);
  origin(-L0/2, -L0/2, );
  init_grid (1 << MAXLEVEL);
  p[inletflow] = dirichlet(1e-3);
  run();
}


event init (i = 0) {
    mask(INLET ? inletflow : none);
    scalar inl[];
    int it = 0;
    do {
        it++;
        foreach()
        {
            f[] = BOOL_CYLINDER ;//* (1. - 1e-3 * (1. + noise()));
            p[] = 1e-3*f[];
            inl[] = INLET;
        }
        if (it>=10) printf("WARNING: does not converge... ");
    }while (adapt_wavelet({f, inl }, (double []){0.01, 0.001},
                          maxlevel = MAXLEVEL, minlevel = 3).nf != 0 && it <= MAXLEVEL+2);

//    p[inletflow] = dirichlet(1e-3);
//    f[inletflow] = dirichlet(1.);
//    boundary (all);
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
  boundary ({f});
  double ff;
  foreach_face() {
    ff = (f[] + f[-1])/2.;
    beta.x[] = -(k / (mu1 + clamp(ff,0,1)*(mu2 - mu1)));
  }

    /**
  inside the circle we add a source via a non-zero divergence.
     */
    foreach() {
        if (BOOL_CYLINDER) {
            zeta[] = exp(-500.*(sq(x) + sq(y)));
            f[] = 1;
        } else
            zeta[] = 0;
    }
    boundary ((scalar *){beta, zeta});
}
#if TREE
event adapt (i++; t <= 50) {
  adapt_wavelet ({f}, (double[]){0.01}, MAXLEVEL);
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

  output_vtu_bin_foreach((scalar *) {l, f, p}, (vector *) {u, beta}, 64, fp, false);
  fclose(fp);
//  @if _MPI
//  if (pid() == 0) {
    sprintf(name, "hs_%4.4d.pvtu", nf);
    sprintf(subname, "hs_%4.4d", nf);
    fp = fopen(name, "w");
    output_pvtu_bin((scalar *) {l, f, p}, (vector *) {u, beta}, 64, fp, subname);
    fclose(fp);
//  }
//  MPI_Barrier(MPI_COMM_WORLD);
//  @endif
    fprintf (ferr, "iteration: %d\n", iteration); fflush (ferr);
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