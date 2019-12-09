/**
# Stokes flow past a periodic array of cylinders
We compare the numerical results with the solution given by the
multipole expansion of [Sangani and Acrivos, 1982](#sangani1982). */
//const vector zerocf[] = {0.,0.,0.};
#define BRINKMAN_PENALIZATION 4
#define DEBUG_BRINKMAN_PENALIZATION 1

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

int maxlevel = 10, minlevel = 4, nc;
double radius;

void calc_solid(scalar fs, vector n_sol, vector target_U){
  vertex scalar phi[];
  face vector face_fs[];
  foreach_vertex() {
    phi[] = (sq(x) + sq(y) < sq(radius) ) ? 1 : -1;
  }
  boundary ({phi});
  fractions (phi, fs, face_fs);
  foreach() {
    foreach_dimension() target_U.x[] = 0;
    n_sol.x[] = x/sqrt(sq(x) + sq(y));
    n_sol.y[] = y/sqrt(sq(x) + sq(y));
  }
  boundary ({fs, target_U, n_sol});
}

int main(int argc, char * argv[]){
    eta_s =1e-5;
    if (argc > 1) {
      eta_s = atof(argv[1]); //convert from string to float
    }
    if (argc > 2) {
      maxlevel = atoi(argv[2]); //convert from string to float
    }
    fprintf(fout, "eta_s=%g maxlevel=%d", eta_s, maxlevel);
    /**
    The domain is the periodic unit square, centered on the origin. */
    L0=1;
    origin (-L0/2., -L0/2.);
    periodic (right);
    periodic (top);

    /**
    We turn off the advection term. The choice of the maximum timestep
    and of the tolerance on the Poisson and viscous solves is not
    trivial. This was adjusted by trial and error to minimize (possibly)
    splitting errors and optimize convergence speed. */

    stokes = true;
    DT = 1e-3;
    TOLERANCE = 1e-8;
    NITERMIN = 5;
    /**
    We do the 9 cases computed by Sangani & Acrivos. The radius is
    computed from the volume fraction. */

    for (nc = 0; nc < 9; nc++) {
        N = 1 << maxlevel;
        radius = sqrt(sq(L0) * sangani[nc][0] / pi);
        run();
    }
}

/**
We need an extra field to track convergence. */

scalar un[];

event init (t = 0){
    int it = 0;
    do {
        calc_solid(fs, n_sol, target_U);
    }while (adapt_wavelet({fs}, (double []){0.001},
                          maxlevel = maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
    /**
    And set acceleration and viscosity to unity. */

    const face vector g[] = {1.,0.};
    a = g;
    mu = fm;

    /**
    We initialize the reference velocity. */
    foreach() un[] = u.y[];
    event("vtk_file");
}

/**
We check for a stationary solution. */

//event logfile (i++; t <= 0.27){
event logfile (i++; i <= 5000){
    double avg = normf_weugene(u.x, fs).avg;
    double du = change_weugene (u.x, un, fs)/(avg + SEPS); //change 1) Linf  2) un = u
    fprintf (fout, "%d %d %d %d %d %d %d %d %.3g %.3g %.3g %.3g %.3g\n",
    maxlevel, i,
    mgp.i, mgp.nrelax, mgp.minlevel,
    mgu.i, mgu.nrelax, mgu.minlevel,
    du, mgp.resa*dt, mgu.resa, statsf_weugene(u.x, fs).sum, normf_weugene(p, fs).max);
    fflush(fout);

//    if(i>1){
//        stats s = statsf_weugene(u.x, fs);
//        double Phi = 1. - s.volume/sq(L0);
//        double Phia = pi*sq(radius)/sq(L0);
//        double U = s.sum/s.volume;
//        double F = sq(L0)/(1. - Phi);
//        fprintf (ferr, "%d %g %g %g %g %d| %g=%g? %g %g\n", maxlevel, sangani[nc][0], F/U, sangani[nc][1], fabs(F/U - sangani[nc][1])/sangani[nc][1], i, Phi, Phia, U, F);
//    }
//    if (t >= 0.267) {
    if (((i > 1) && (du < 1e-3))|| (i == 5000)) {
        /**
        We output the non-dimensional force per unit length on the
        cylinder $F/(\mu U)$, together with the corresponding value from
        Sangani & Acrivos and the relative error. */
        stats s = statsf_weugene(u.x, fs);
        double Phi = 1. - s.volume/sq(L0);
        double Phia = pi*sq(radius)/sq(L0);
        double U = s.sum/s.volume;
        double F = sq(L0)/(1. - Phi);
        fprintf (ferr, "%d %g %g %g %g %d| %g=%g? F:%g dt:%g t:%g U:%g\n", maxlevel, sangani[nc][0], F/U, sangani[nc][1], fabs(F/U - sangani[nc][1])/sangani[nc][1], i, Phi, Phia, F, dt, t, U);
//        stats s = statsf(u);
//        double Phi = 1. - s.volume/sq(L0);
//        double U = s.sum/s.volume;
//        double F = sq(L0)/(1. - Phi);
//        fprintf (ferr,
//        "%d %g %g %g %g %d\n", maxlevel, sangani[nc][0], F/U, sangani[nc][1],
//        fabs(F/U - sangani[nc][1])/sangani[nc][1], i);

//        view (fov = 9.78488, tx = 0.250594, ty = -0.250165);
//        draw_vof ("fs",  lc = {1,0,0}, lw = 2); // draw line lc -color, lw -width
//        squares ("u.x", linear = 1, spread = -1); // spread<0 => color is distributed min max
//        cells();
//        char subname[80]; sprintf(subname, "mesh-%d.png", nc);
//        save (subname);
//
        fprintf(fout, "stationary flow nc = %d i = %d du = %g", nc, i, du);
        fflush(fout);
        return 9;
    }
}

//Output
#include "../src_local/output_vtu_foreach.h"
event vtk_file (t += 0.01){
    char subname[80]; sprintf(subname, "br");
    scalar l[];
    vorticity (u, omega);
    foreach() {l[] = level; omega[] *= 1 - fs[];}
    output_vtu_MPI( (scalar *) {l, omega, fs, p}, (vector *) {u, uf, dbp, utau, grad_utau_n, n_sol, target_U}, subname, L0/pow(2, minlevel));
}

#define ADAPT_SCALARS {fs, omega}
#define ADAPT_EPS_SCALARS {1e-3, 1e-2}
event adapt (i++){
    double eps_arr[] = ADAPT_EPS_SCALARS;
//    MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
    calc_solid(fs, n_sol, target_U);
}

/**
The non-dimensional drag force per unit length closely matches the
results of Sangani & Acrivos. For $\Phi=0.75$ and level 8 there is
only about 6 grid points in the width of the gap between cylinders.

~~~gnuplot Non-dimensional drag force per unit length
set xlabel 'Volume fraction'  font ",10"
set ylabel 'k_0'  font ",10"
set logscale y
set grid
set key top left
set tics font "Helvetica,10"
set format y "10^{%L}"
plot '< grep "^10" log_1e-6' u 2:4 ps 5 lw 5 t 'Sangani and Acrivos, 1982','' u 2:3 ps 5 pt 6 lw 5 t '10 levels'
~~~
 plot '< grep "^11" log_1e-6'  u 2:3 w lp ps 1 pt 6 lw 2 t 'BP eta=1e-6, 11 levels','< grep "^11" log_1e-5' u 2:3 w lp ps 1 pt 30 lw 2 t 'BP eta=1e-5, 11 levels','< grep "^11" log_1e-4' u 2:3 w lp ps 1 pt 33 lw 2 t 'BP eta=1e-4, 11 levels','< grep "^11" log_1e-3' u 2:3 w lp ps 1 pt 10 lw 2 t 'BP eta=1e-3, 11 levels','< grep "^11" log_1e-2' u 2:3 w lp ps 1 pt 12 lw 2 t 'BP eta=1e-2, 11 levels','< grep "^11" log_1e-6' u 2:4 w lp pt 2 ps 1 lc rgb "orange-red" lw 2 t 'Sangani and Acrivos, 1982'
 
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