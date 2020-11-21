/**
# Stokes flow past a periodic array of cylinders
We compare the numerical results with the solution given by the
multipole expansion of [Sangani and Acrivos, 1982](#sangani1982). */
#define BRINKMAN_PENALIZATION 1
#define DEBUG_BRINKMAN_PENALIZATION 1
#define RELATIVE_RESIDUAL
#undef SEPS
#define SEPS 1e-30
#define DEBUG_MINMAXVALUES
#define DEBUG_OUTPUT_VTU_MPI

scalar fs[], omega[];
vector Us[];
scalar divutmp[], deltap[];
vector conv_term[], deltag[];
const face vector grav[] = {1.,0.};
#include "../src_local/centered-weugene.h"
#include "view.h"
#include "../src_local/output_vtu_foreach.h"

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

int maxlevel = 8, minlevel = 4, nc;
double radius, epsthresh=1e-2;

void calc_solid(scalar fs){
  vertex scalar phi[];
  face vector face_fs[];
  foreach_vertex() {
    phi[] = sq(radius) - sq(x) - sq(y);
  }
  boundary ({phi});
  fractions (phi, fs, face_fs);
  boundary ({fs});
}

int main(int argc, char * argv[]){
    eta_s =1e-15;
    if (argc > 1) {
      eta_s = atof(argv[1]); //convert from string to float
    }
    if (argc > 2) {
      maxlevel = atoi(argv[2]); //convert from string to float
    }
	if (argc > 3) {
		epsthresh = atof(argv[3]); //convert from string to float
	}
    fprintf(fout, "eta_s=%g maxlevel=%d\n", eta_s, maxlevel);
    /**
    The domain is the periodic unit square, centered on the origin. */
    L0 = 1;
    origin (-L0/2., -L0/2.);
    periodic (right);
    periodic (top);
    /**
    And set acceleration and viscosity to unity. */

    a = grav;
    mu = fm;
    /**
    We turn off the advection term. The choice of the maximum timestep
    and of the tolerance on the Poisson and viscous solves is not
    trivial. This was adjusted by trial and error to minimize (possibly)
    splitting errors and optimize convergence speed. */

    stokes = true;
    DT = 1e-3;
    TOLERANCE = 1e-8;
//    NITERMIN = 1;
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
        calc_solid(fs);
    }while (adapt_wavelet({fs}, (double []){1e-3},
                          maxlevel = maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
    /**
    We initialize the reference velocity. */
    foreach() {
        un[] = u.y[];
        g.x[] = grav.x[];
    }
    //event("vtk_file");
}

event end_timestep(i++)
{
    double eps_arr[] = {1, 1, 1, 1, 1, 1, 1};
    MinMaxValues((scalar *){p, u.x, u.y, g.x, g.y}, eps_arr);
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
    if (i>1){
        /**
            We output the non-dimensional force per unit length on the
            cylinder $F/(\mu U)$, together with the corresponding value from
            Sangani & Acrivos and the relative error. */
        stats s = statsf_weugene(u.x, fs);
        double Phi = 1. - s.volume / sq(L0);
        double Phia = pi * sq(radius) / sq(L0);
        double U = s.sum / s.volume;
        double F = sq(L0) / (1. - Phi);
        fprintf (ferr,
        "+++%d vol.fr.theory=%g F/U=%g F/Uth=%g delta(F/U-theory)=%g i=%d ifp=%d| %g=%g? F:%g dt:%g t:%g U:%g\n", maxlevel, sangani[nc][0], F/U, sangani[nc][1],
        fabs(F/U - sangani[nc][1])/sangani[nc][1], i, iter_fp, Phi, Phia, F, dt, t, U);
        fprintf (fout, "stationary flow nc = %d i = %d du = %g", nc, i, du);
        fflush (fout);

        if (((i > 1) && (du/dt < 1e-2))|| (i == 5000)) {
            fprintf (ferr,
            "%d vol.fr.theory=%g F/U=%g F/Uth=%g delta(F/U-theory)=%g i=%d ifp=%d| %g=%g? F:%g dt:%g t:%g U:%g\n", maxlevel, sangani[nc][0], F/U, sangani[nc][1],
            fabs(F/U - sangani[nc][1])/sangani[nc][1], i, iter_fp, Phi, Phia, F, dt, t, U);
            event("vtk_file");
            return 9;
        }
    }

}

//Output
event vtk_file (t += 0.1){
    char subname[80]; sprintf(subname, "cyl_br");
    scalar l[];
    vorticity (u, omega);
    foreach() {l[] = level; omega[] *= 1 - fs[];}
//    output_vtu_MPI( (scalar *) {l, omega, fs, p}, (vector *) {u, uf, dbp, utau, grad_utau_n, n_sol, target_U}, subname, L0/pow(2, minlevel));
    output_vtu_MPI( (scalar *) {l, omega, fs, p}, (vector *) {u, dbp, utau, total_rhs, conv_term, g},(vector *) {uf}, subname, t+dt);
}

#define ADAPT_SCALARS {u.x, u.y}
#define ADAPT_EPS_SCALARS {epsthresh, epsthresh}
event adapt (i++){
    double eps_arr[] = ADAPT_EPS_SCALARS;
    MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
	calc_solid(fs);
    fs.refine = fs.prolongation = fraction_refine;
    boundary({fs});
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
plot '< grep "^10" log' u 2:4 ps 5 lw 5 t 'Sangani and Acrivos, 1982','' u 2:3 ps 5 pt 6 lw 5 t '10 levels'
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
