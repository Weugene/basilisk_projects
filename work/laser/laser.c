/**
# Wave equation is
$$
i\frac{\psi}{t}=-\frac{\nabla^2\psi}{2M} + \frac12 kx^2\psi
$$

Wave function $\psi$ is $a+ib$, after splitting into real and imaginary parts we have:
$$
\partial_t a = -\frac{\nabla^2 b}{2M} + \frac12 kx^2 b
$$
$$
\partial_t b = \frac{\nabla^2 a}{2M} - \frac12 kx^2 a
$$

We will use a Cartesian (multi)grid, the generic time loop and the
time-implicit diffusion solver. */

#include "grid/multigrid.h"
#include "run.h"
#include "diffusion.h"

/**
We need scalar fields for the imaginary and real part of wave function. */

scalar a[], b[];


double Mass = 1., ka = 4.5, A0 = 1.;

/**
The generic time loop needs a timestep. We will store the statistics
on the diffusion solvers in `mgd1` and `mgd2`. */

double dt;
mgstats mgd1, mgd2;

/**
## Parameters

We change the size of the domain `L0` and set the tolerance of the
implicit diffusion solver. */

int main()
{
  N=1024;
  init_grid (N);
  size (N);
  origin(-L0/2, -L0/2, -L0/2);
  TOLERANCE = 1e-4;

  run();
}

/**
## Initial conditions */

event init (i = 0)
{

  foreach() {
    double r2 = 0.0;
    foreach_dimension(){r2 += x*x;}
    a[] = A0 * exp(-r2 / (2.0*w0*w0));
    b[] = 0.0;
  }
  boundary ({a, b});
}

/**
## Outputs

Here we create an mpeg animation of the $C_1$ concentration. The
`spread` parameter sets the color scale to $\pm$ twice the standard
deviation. */

event movie (i = 1; i += 10)
{
  scalar modulus[];
  foreach(){
    modulus[] = sqrt(a[]*a[] + b[]*b[]);
  }
  output_ppm (modulus, linear = true, spread = 2, file = "f.ppm_");
  fprintf (stderr, "%d %g %g %d %d\n", i, t, dt, mgd1.i, mgd2.i);
}

/**
We make a PNG image of the final "pseudo-stationary" solution. */

event final (t = 3000)
{
  char name[80];
  output_ppm (a, file = name, n = N, linear = true, spread = 2);
}

/**
## Time integration */

event integration (i++)
{
  /**
  We first set the timestep according to the timing of upcoming
  events. We choose a maximum timestep of 1 which ensures the stability
  for this example. */

  dt = dtnext (1.);
/**
  $$
  \partial_t a = -\frac{\nabla^2 b}{2M} + \frac12 kx^2 b
  $$
  $$
  \partial_t b = \frac{\nabla^2 a}{2M} - \frac12 kx^2 a
$$*/
  scalar r[];

  foreach() {
    double r2 = 0;
    foreach_dimension(){r2 += x*x;}
    r[] = 0.5*k*r2*b;
  }
  mgd1 = diffusion (a, dt, r = r);
  foreach() {
    r[] = k*kb*a[];
  }
  const face vector c[] = {0, 0};
  mgd2 = diffusion (b, dt, c, r);
}

/**
## Results

We get the following stable [Turing
patterns](http://en.wikipedia.org/wiki/The_Chemical_Basis_of_Morphogenesis).

|:------:|:-----:|:-----:|
| ![](brusselator/mu-0.04.png) | ![](brusselator/mu-0.1.png) | ![](brusselator/mu-0.98.png) |
| $\mu=0.04$                   |    $\mu=0.1$ (stripes)      |    $\mu=0.98$ (hexagons)     |

  : [Animation of the transitions](brusselator/f.mp4)
*/
