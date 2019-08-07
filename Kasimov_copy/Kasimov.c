
/**
# Atomisation of a pulsed liquid jet

A dense cylindrical liquid jet is injected into a stagnant lighter
phase (density ratio 1/27.84). The inflow velocity is modulated
sinusoidally to promote the growth of primary shear
instabilities. Surface tension is included and ultimately controls the
characteristic scale of the smallest droplets.

We solve the two-phase Navier--Stokes equations with surface
tension. We need the *tag()* function to count the number of
droplets. We generate animations online using Basilisk View. */

#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "vtk.h"
#include "view.h"
/**
We define the radius of the jet, the initial jet length, the Reynolds
number and the surface tension coefficient. */

#define LEVEL 10
#define length 0.1
#define feps 0.001
#define ueps 0.001
#define RE 5000.0
#define CA 1.0e-2
//#define FROUDE 0.4
#define RHOR 100.
#define MUR 100.
/**
The default maximum level of refinement is 10 and the error threshold
on velocity is 0.1. */

int MAXLEVEL = 12;

static int iteration=0;
/**
To impose boundary conditions on a disk we use an auxilliary volume
fraction field *f0* which is one inside the cylinder and zero
outside. We then set an oscillating inflow velocity on the
left-hand-side and free outflow on the right-hand-side. */

scalar f0[];
u.n[left]  = dirichlet(f0[]);
//u.n[left]  = dirichlet(f0[]*(1. + 0.1*sin (10.*2.*pi*t)));
u.t[left]  = dirichlet(0);
#if dimension > 2
u.r[left]  = dirichlet(0);
#endif
p[left]    = neumann(0);
f[left]    = neumann(0.0);

u.n[right] = neumann(0);
//u.t[right] = neumann(0);
p[right]   = dirichlet(0);
f[right]   = neumann(0);


/**
The program can take two optional command-line arguments: the maximum
level and the error threshold on velocity. */

int main (int argc, char * argv[])
{
    if (argc > 1)
        MAXLEVEL = atoi (argv[1]);
//    if (argc > 2)
//        ueps = atof (argv[2]);

    /**
    The initial domain is discretised with $64^3$ grid points. We set
    the origin and domain size. */
    init_grid (1 <<LEVEL);
    origin (0, 0, 0);
    size (1.);
    periodic (top);
    /**
    We set the density and viscosity of each phase as well as the
    surface tension coefficient and start the simulation. */

    rho1 = 1.; rho2 = 1./RHOR;
    mu1 = 1./RE; mu2 = 1./(MUR*RE);
    f.sigma = 1./(RE*CA);
    //G.y = - 1./sq(FROUDE);
    run();
}

/**
## Initial conditions */

double front (double x, double y, double z)
{
    double phi = HUGE;
    return min(phi, 0.2*0.5*(1.0 - tanh((fabs(y-L0/2.0) - 0.02)/0.01)) +0.0-x);

    //return min(phi, -fabs(y-0.5)+(length - x));
    //return min(phi, -0.005*sin(10.0*2.0*pi*y)+(length - x));
}


event init (t = 0) {

    if (!restore (file = "restart")) {

        /**
        We use a static refinement down to *maxlevel* in a cylinder 1.2
        times longer than the initial jet and twice the radius. */
        int it = 0;
        do {
            it++;
            fraction(f0, front (x, y, z));
            if (it>=10) printf("WARNING: does not converge... ");
        }while (adapt_wavelet({f0, u.x}, (double []){feps, ueps},
                              maxlevel = MAXLEVEL, minlevel = 5).nf != 0 && it <= 10);


        /**
        We then use this to define the initial jet and its velocity. */

        foreach() {
            f[] = f0[];
            u.x[] = f0[];
        }
        boundary ({f,u.x});
    }
}

/**
## Outputs

We log some statistics on the solver. */

event logfile (i++) {
    if (i == 0)
        fprintf (ferr,
                 "t dt mgp.i mgpf.i mgu.i grid->tn perf.t perf.speed\n");
    fprintf (ferr, "%g %g %d %d %d %ld %g %g\n",
             t, dt, mgp.i, mgpf.i, mgu.i,
             grid->tn, perf.t, perf.speed);
}

/**
We generate an animation using Basilisk View. */

event movie (t += 1e-1; t <= 10)
{
#if dimension == 2
  scalar omega[];
    vorticity (u, omega);
  view (tx = 0.5);
  clear();
  draw_vof ("f");
  squares ("omega", linear = true, spread = 10);
  box ();
#else // 3D
    scalar pid[];
    foreach()
    pid[] = fmod(pid()*(npe() + 37), npe());
    boundary ({pid}); // not used for the moment
    view (camera = "iso",
          fov = 14.5, tx = -0.418, ty = 0.288,
          width = 1600, height = 1200);
    clear();
    draw_vof ("f");
#endif // 3D
    save ("movie.mp4");
}
//Don't run in parallel!!! it does not work in parallel
#if _MPI != 1
event vtk_file (i += 100)
{
    char name[80];
    sprintf(name, "list.vtk.%d", iteration++);
    printf("name %s", name);
    FILE * fplist = fopen (name, "w");
    scalar l[];
    foreach()
        l[] = level;
    output_vtk ({f, u.x, u.y, rhov, p, l}, 512, fplist, false );
}
#endif
/**
We save snapshots of the simulation at regular intervals to
restart or to post-process with [bview](/src/bview). */

event snapshot (i += 500) {
    char name[80];
    sprintf (name, "dump-%d", i);
    scalar pid[];
    foreach()
        pid[] = fmod(pid()*(npe() + 37), npe());
    boundary ({pid});
    dump (name);
}



/**
## Mesh adaptation

We adapt the mesh according to the error on the volume fraction field
and the velocity. */

event adapt (i++) {
    adapt_wavelet ({f, u}, (double[]){feps,ueps,ueps,ueps}, MAXLEVEL);
}

/**
## Running in parallel

### On Occigen

To run in 3D on
[occigen](https://www.cines.fr/calcul/materiels/occigen/), we can do on
the local machine

~~~bash
local% qcc -source -grid=octree -D_MPI=1 atomisation.c
local% scp _atomisation.c occigen.cines.fr:
~~~

and on occigen (to run on 64*24 = 1536 cores, with 12 levels of refinement)

~~~bash
occigen% sbatch --nodes=64 --time=5:00:00 run.sh
~~~

with the following `run.sh` script

~~~bash
#!/bin/bash
#SBATCH -J basilisk
#SBATCH --nodes=1
#SBATCH --constraint=HSW24
#SBATCH --ntasks-per-node=24
#SBATCH --threads-per-core=1
#SBATCH --time=10:00
#SBATCH --output basilisk.output
#SBATCH --exclusive

LEVEL=12

module purge
module load openmpi
module load intel

NAME=atomisation
mpicc -Wall -std=c99 -O2 _$NAME.c -o $NAME \
      -L$HOME/gl -lglutils -lfb_osmesa -lOSMesa -lGLU -lm
srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS ./$NAME $LEVEL \
     2> log-$LEVEL-$SLURM_NTASKS > out-$LEVEL-$SLURM_NTASKS
~~~

Note that this assumes that the gl libraries have been
[installed](/src/gl/INSTALL) in `$HOME/gl`.

### On Mesu

To run in 3D on [mesu](http://iscd.upmc.fr/expertise/mesu), we can do
on the local machine

~~~bash
local% qcc -source -grid=octree -D_MPI=1 atomisation.c
local% scp _atomisation.c mesu.dsi.upmc.fr:
~~~

and on mesu (to run on 672 cores, with 12 levels of refinement)

~~~bash
mesu% qstat -u popinet
mesu% qsub run.sh
~~~

with the following `run.sh` script

~~~bash
#!/bin/bash
#PBS -l select=28:ncpus=24:mpiprocs=24
#PBS -l walltime=12:00:00
#PBS -N atomisation
#PBS -j oe
# load modules
module load mpt
mpicc -Wall -O2 -std=c99 _atomisation.c -o atomisation \
     -L$HOME/gl -lglutils -lfb_osmesa -lOSMesa -lGLU -lm
# change to the directory where program job_script_file is located
cd $PBS_O_WORKDIR
# mpirun -np 672 !!!! does not work !!!!
mpiexec_mpt -n 672 ./atomisation 12 2>> log >> out
~~~

Note that this assumes that the gl libraries have been
[installed](/src/gl/INSTALL) in `$HOME/gl`.

## Results

![Atomisation of a pulsed liquid jet. 4096^3^ equivalent
 resolution.](atomisation/movie.mp4)(width="800" height="600")

The gain in number of grid points of the adaptive simulation (relative
to a constant 4096^3^ Cartesian grid discretisation) is illustrated
below as well as the computational speed (in
points.timesteps/sec). The simulation on Mesu is restarted at $t=2.6$
(after 12 hours of runtime).

~~~gnuplot Compression ratio and computational speed
set xlabel 'Time'
set logscale y
set key center right
plot 'atomisation.occigen' u 1:8 w l t 'speed (occigen 1536 cores)', \
     'atomisation.mesu' u 1:8 w l t 'speed (mesu 672 cores)', \
     'atomisation.mesu' u 1:(4096.**3./$6) w l t 'compression ratio'
~~~

~~~gnuplot Histogram of droplet volumes at $t=3.7$
reset
set xlabel 'log10(Volume)'
set ylabel 'Number of droplets'
binstart = -100
binwidth = 0.2
set boxwidth binwidth
set style fill solid 0.5
Delta = 3./2**12 # minimum cell size
minvol = log10(Delta**3) # minimum cell volume
set arrow from minvol, graph 0 to minvol, graph 1 nohead
set label "Minimum cell volume" at minvol - 0.2, graph 0.35 rotate by 90
t = "3.7"
plot "< awk '{if ($2==".t.") print $0;}' atomisation.out" u \
     (binwidth*(floor((log10($4)-binstart)/binwidth)+0.5)+binstart):(1.0) \
     smooth freq w boxes t ''
~~~

## See also

* [Same example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/atomisation.html)
*/
