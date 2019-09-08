/**
# Convergence of the Runge--Kutta solvers

We solve numerically
$$
\frac{du}{dt} = tu
$$
with the initial condition $u(0)=1$. The solution is $u(t)=e^{t^2/2}$.
*/
//#include "predictor-corrector.h"
#include "runge-kutta.h"
#include "run.h"

/**
This function returns the right-hand-side of the equation i.e. $tu$. */

static void du (scalar * ul, double t, scalar * kl)
{
  scalar u = ul[0], k = kl[0];
  foreach()
    k[] = t*u[];
}
scalar u[], uerr[];
vector uvec[];
int order=1;


double emax = 0.;

int main()
{
  init_grid (256);
  order = 1;
  dt = 1e-2;
  run();
  printf ("\n");
  fprintf (stderr, "%g %g %d\n", dt, emax, order);
}
event init(t = 0){
    foreach()
    u[] = 1.; // the initial condition
    emax = 0.;
    foreach()
    foreach_dimension()
    uvec.x[] =1;
}
event stepik(i++){
    foreach() {
        uerr[] = fabs (u[] - exp(t*t/2.));
        if (uerr[] > emax)
            emax = uerr[];
        printf ("%g %g %g\n", t, u[], uerr[]);
    }
    runge_kutta ({u}, t, dt, du, order);

}
//Output
static int iteration=0;
#include "output_fields/output_vtu_foreach.h"
event vtk_file (i++;t<1)
{
//    int nf = iteration;
//    scalar l[];
//    foreach()
//        l[] = level;
//
//    char name[80], subname[80];
//    FILE *fp;
//    sprintf(name, "hs%4.4d_n%3.3d.vtu", nf, pid());
//    fp = fopen(name, "w");
//
//    output_vtu_bin_foreach((scalar *) {u, uerr, l}, (vector *) {}, 64, fp, false);
//
//    fclose(fp);
//    @if _MPI
//    if (pid() == 0) {
//        sprintf(name, "hs%4.4d.pvtu", nf);
//        sprintf(subname, "hs%4.4d", nf);
//        fp = fopen(name, "w");
//        output_pvtu_bin((scalar *) {u, uerr, l}, (vector *) {}, 64, fp, subname);
//        fclose(fp);
//    }
//    MPI_Barrier(MPI_COMM_WORLD);
//    @endif
//    fprintf (ferr, "iteration: %d\n", iteration); fflush (ferr);
//    iteration++;
}

/**
~~~gnuplot Error convergence for different orders
set xlabel 'dt'
set ylabel 'error'
set logscale
set key top left
set ytics format '%.0e'
plot "< grep '1$' log" pt 7 t '', 15.*x t '15 dt',  \
     "< grep '2$' log" pt 7 t '', 4.*x*x t '4 dt^2', \
     "< grep '4$' log" pt 7 t '', x**4/2. t 'dt^4/2'
~~~
*/
