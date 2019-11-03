/**
# Stokes flow past a periodic array of cylinders

We compare the numerical results with the solution given by the
multipole expansion of [Sangani and Acrivos, 1982](#sangani1982). */

#define BRINKMAN_PENALIZATION
#define DEBUG_BRINKMAN_PENALIZATION

#undef SEPS
#define SEPS 1e-30
//#define DEBUG_MINMAXVALUES
//#define DEBUG_OUTPUT_VTU_MPI

scalar fs[], omega[];
vector Us[];

#include "../src_local/centered-weugene.h"
#include "view.h"

/**
This is Table 1 of [Sangani and Acrivos, 1982](#sangani1982), where
the first column is the volume fraction $\Phi$ of the cylinders and
the second column is the non-dimensional drag force per unit length of
cylinder $F/(\mu U)$ with $\mu$ the dynamic vicosity and $U$ the
average fluid velocity. */

static double sangani[9][2] = {
        {0.05, 15.56},
        {0.10, 24.83},
        {0.20, 51.53},
        {0.30, 102.90},
        {0.40, 217.89},
        {0.50, 532.55},
        {0.60, 1.763e3},
        {0.70, 1.352e4},
        {0.75, 1.263e5}
};

/**
We will vary the maximum level of refinement, *nc* is the index of the
case in the table above, the radius of the cylinder will be computed
using the volume fraction $\Phi$. */

int maxlevel = 10, minlevel =4, nc;
double radius;

int main(int argc, char * argv[]){
    eta_s =1e-5;
    if (argc > 1) {
      eta_s = atof(argv[1]); //convert from string to float
    }
    /**
    The domain is the periodic unit square, centered on the origin. */

    size (1.);
    origin (-L0/2., -L0/2.);
    periodic (right);
    periodic (top);

    /**
    We turn off the advection term. The choice of the maximum timestep
    and of the tolerance on the Poisson and viscous solves is not
    trivial. This was adjusted by trial and error to minimize (possibly)
    splitting errors and optimize convergence speed. */

    stokes = true;
    DT = 1e-4;
    TOLERANCE = 1e-5;
    NITERMIN = 5;
    /**
    We do the 9 cases computed by Sangani & Acrivos. The radius is
    computed from the volume fraction. */

    for (nc = 0; nc < 9; nc++) {
        N = 1 << 8;
        radius = sqrt(sq(L0) * sangani[nc][0] / pi);
        run();
    }
}

/**
We need an extra field to track convergence. */

scalar un[];

event init (t = 0)
{
    /**
    We initialize the embedded geometry. */
    foreach() fs[] = sq(x) + sq(y) < sq(radius);
    boundary ({fs});
    /**
    And set acceleration and viscosity to unity. */

    const face vector g[] = {1.,0.};
    a = g;
    mu = fm;
    /**
    We initialize the reference velocity. */
    foreach()
    un[] = u.y[];
}

/**
We check for a stationary solution. */

event logfile (i++; i <= 5000){
  scalar u_tmp[];
  foreach() u_tmp[] = u.x[]*(1-fs[]);
  double avg = normf(u_tmp).avg;
  double du = change (u_tmp, un)/(avg + SEPS); //change 1) Linf  2) un = u
//  fprintf (fout, "%d %d %d %d %d %d %d %d %.3g %.3g %.3g %.3g %.3g\n",
//  maxlevel, i,
//  mgp.i, mgp.nrelax, mgp.minlevel,
//  mgu.i, mgu.nrelax, mgu.minlevel,
//  du, mgp.resa*dt, mgu.resa, statsf(u.x).sum, normf(p).max);
  if (i > 1 && du < 1e-3 || i == 5000) {
  /**
  We output the non-dimensional force per unit length on the
  cylinder $F/(\mu U)$, together with the corresponding value from
  Sangani & Acrivos and the relative error. */

  stats s = statsf(u_tmp);
  double Phi = 1. - s.volume/sq(L0);
  double U = s.sum/s.volume;
  double F = sq(L0)/(1. - Phi);
  fprintf (ferr,
  "%d %g %g %g %g\n", maxlevel, sangani[nc][0], F/U, sangani[nc][1],
  fabs(F/U - sangani[nc][1])/sangani[nc][1]);

  view (fov = 9.78488, tx = 0.250594, ty = -0.250165);
  draw_vof ("fs",  lc = {1,0,0}, lw = 2); // draw line lc -color, lw -width
  squares ("u.x", linear = 1, spread = -1); // spread<0 => color is distributed min max
  cells();
  char subname[80]; sprintf(subname, "mesh-%d.png", nc);
  save (subname);

  fprintf(fout, "stationary flow nc = %d i = %d du = %g", nc, i, du);
  return 9;
  }
}

//#define ADAPT_SCALARS {fs, u}
//#define ADAPT_EPS_SCALARS {1e-3,2e-5,2e-5}

#define ADAPT_SCALARS {fs, omega}
#define ADAPT_EPS_SCALARS {1e-3,1e-3}
event adapt (i++){
  double eps_arr[] = ADAPT_EPS_SCALARS;
//  MinMaxValues(ADAPT_SCALARS, eps_arr);
  adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
  foreach() fs[] = sq(x) + sq(y) < sq(radius);
}
//Output
#include "../src_local/output_vtu_foreach.h"
event end_timestep (t += 0.01){
  char subname[80]; sprintf(subname, "br");
  scalar l[], omega[];
  vorticity (u, omega);
  foreach() l[] = level;
  output_vtu_MPI( (scalar *) {l, omega, fs, p}, (vector *) {u, dbp}, subname, L0/pow(2, minlevel));
}
/**
The non-dimensional drag force per unit length closely matches the
results of Sangani & Acrivos. For $\Phi=0.75$ and level 8 there is
only about 6 grid points in the width of the gap between cylinders.

~~~gnuplot Non-dimensional drag force per unit length
set xlabel 'Volume fraction'
set ylabel 'k_0'
set logscale y
set grid
set key top left
plot '< grep "^8" log' u 2:4 ps 1 lw 2 t 'Sangani and Acrivos, 1982', \
     '' u 2:3 ps 1 pt 6 lw 2 t '8 levels'
~~~

This can be further quantified by plotting the relative error. It
seems that at low volume fractions, the error is independent from the
mesh refinement. This may be due to other sources of errors, such as
the splitting error in the projection scheme. This needs to be
explored further.

~~~gnuplot Relative error on the drag force
set ylabel 'Relative error'
plot '< grep "^6" log' u 2:5 w lp t '6 levels',  \
     '< grep "^7" log' u 2:5 w lp t '7 levels', \
     '< grep "^8" log' u 2:5 w lp t '8 levels', \
     '< grep "^9" log' u 2:5 w lp t '9 levels'
~~~

The adaptive mesh for 9 levels of refinement illustrates the automatic
refinement of the narrow gap between cylinders.

![Adaptive mesh at level 9 for $\Phi=0.75$.](cylinders/mesh.png)

## References

~~~bib
@article{sangani1982,
  title={Slow flow past periodic arrays of cylinders 
  with application to heat transfer},
  author={Sangani, AS and Acrivos, A},
  journal={International Journal of Multiphase Flow},
  volume={8},
  number={3},
  pages={193--206},
  year={1982},
  publisher={Elsevier}
}
~~~
*/