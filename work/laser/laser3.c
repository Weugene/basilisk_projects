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
// #include "grid/quadtree.h"
/**
We need scalar fields for the imaginary and real part of wave function. */

scalar a[], b[], modulus[];
// FILE * fp;

double Mass = 1., ka = 10, A0 = 1., w0=.1;

/**
The generic time loop needs a timestep. We will store the statistics
on the diffusion solvers in `mgd1` and `mgd2`. */

double dt;
mgstats mgd1, mgd2;

/**
## Parameters

We change the size of the domain `L0` and set the tolerance of the
implicit diffusion solver. */
// a[right] = dirichlet(0);
// b[right] = dirichlet(0);
// a[left] = neumann(0);
// b[left] = neumann(0);
// a[top] = neumann(0);
// b[top] = neumann(0);
// a[bottom] = neumann(0);
// b[bottom] = neumann(0);
int main()
{

  N=512;
  init_grid (N);
  size (1);
  origin(-L0/2, -L0/2, -L0/2);
  periodic(right);
  periodic(top);
  // fp = popen("ppm2mp4 movie.mp4", "w");
  TOLERANCE = 1e-3;

  run();
}

/**
## Initial conditions */

event init (i = 0)
{
  foreach() {
    double r2 = sq(x) + sq(y);
    a[] = A0 * exp(-r2 / (2.0*w0*w0));
    b[] = A0 * exp(-r2 / (2.0*w0*w0));
  }
  boundary ({a, b});
  event ("vtk_file");
  printf("dim=%d", dimension);
}

/**
## Outputs

Here we create an mpeg animation of the $C_1$ concentration. The
`spread` parameter sets the color scale to $\pm$ twice the standard
deviation. */

// event movie (i++)
// {
  // foreach(){
  //   modulus[] = sqrt(sq(a[]) + sq(b[]));
  // }
  // output_ppm(modulus, fp = fp, n = N);
  // output_ppm (modulus, linear = true, spread = 2, file = "f.ppm");
  // output_ppm (modulus, linear = true, spread = 2, file = "f.ppm_");
// }


/**
## Time integration */
/**
This function returns the right-hand-side of the equation */
/**
  $$
  \partial_t a = -\frac{\nabla^2 b}{2M} + \frac12 kx^2 b
  $$
  $$
  \partial_t b = \frac{\nabla^2 a}{2M} - \frac12 kx^2 a
$$*/

event integration (i++,last)
{
  /**
  We first set the timestep according to the timing of upcoming
  events. We choose a maximum timestep of 1 which ensures the stability
  for this example. */

  dt = dtnext (1.0e-6);
  scalar r[];
  const face vector c[] = {0, 0};
  double lap_a, lap_b, r2;
  foreach() {
    r2 = sq(x) + sq(y);
    lap_b = 0;
    foreach_dimension(){
      lap_b += (b[1] - 2.0*b[] + b[-1])/sq(Delta);
    }
    r[] = 0.5*ka*r2*b[] - (0.5/Mass)*lap_b;
  }
  mgd1 = diffusion (a, dt, c, r);
  boundary ({a});

  foreach() {
    r2 = sq(x) + sq(y);
    lap_a = 0;
    foreach_dimension(){
      lap_a += (a[1] - 2.0*a[] + a[-1])/sq(Delta);
    }
    r[] = -0.5*ka*r2*a[] + (0.5/Mass)*lap_a;
  }
  mgd2 = diffusion (b, dt, c, r);
  boundary ({b});
}

/**
Output*/
int iteration = 0;
// #include "output_fields/output_vtu_foreach.h"
// event vtk_file (i ++; t<=1)
// {
//   int nf = iteration;
//   scalar l[];
//   scalar modulus[];
//   foreach(){
//     // l[] = level;
//     modulus[] = sqrt(sq(a[]) + sq(b[]));
//   }
//
//   char name[80], subname[80];
//   FILE *fp;
//   vector vec[];
//   sprintf(name, "hs_%4.4d_n%3.3d.vtu", nf, pid());
//   fp = fopen(name, "w");
//   output_vtu_bin_foreach((scalar *) {l, modulus, a, b}, (vector *){vec}, 64, fp, false);
//   fclose(fp);
//   sprintf(name, "hs_%4.4d.pvtu", nf);
//   sprintf(subname, "hs_%4.4d", nf);
//   fp = fopen(name, "w");
//   output_pvtu_bin((scalar *) {l, modulus, a, b}, (vector *){vec}, 64, fp, subname);
//   fclose(fp);

// }
// #include "vtk.h"
// event images (t=0.01;t<10) {
//    foreach(){
//      modulus[] = sqrt(sq(a[]) + sq(b[]));
//    }
//    static FILE * fpff = fopen ("modulus.ppm", "w");
//    output_ppm (modulus, fpff, min = 0, max = 1);
//    char name[80];
//    printf("iter=%d\n", iteration);
//    sprintf(name, "list.vtk.%d", iteration);
//    printf("name %s", name);
//    FILE * fplist = fopen (name, "w");
//    output_vtk ({a, b, modulus}, N, fplist,  true);
//    iteration++;
// }

#include "vtk.h"
void vtk_file (){
   // static FILE * fpff = fopen ("modulus.ppm", "w");
   // output_ppm (modulus, fpff, min = 0, max = 1);
   char name[80];
   sprintf(name, "list.vtk.%d", iteration);
   printf("name %s", name);
   FILE * fplist = fopen (name, "w");
   output_vtk ({a, b, modulus}, N, fplist,  true);
   iteration++;
}
