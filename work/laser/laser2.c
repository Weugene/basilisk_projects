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

// #include "grid/multigrid.h"
// #include "run.h"
// #include "diffusion.h"
#include "runge-kutta.h"
/**
We need scalar fields for the imaginary and real part of wave function. */
void vtk_file ();
scalar a[], b[], modulus[];
FILE * fp;

double Mass = 1., ka = 4.5, A0 = 1., w0=.1;
int iteration = 0;
/**
The generic time loop needs a timestep. We will store the statistics
on the diffusion solvers in `mgd1` and `mgd2`. */

double dt;
static void du (scalar * ul, double t, scalar * rhs)
{
  scalar a = ul[0], b = ul[1], a_rhs = rhs[0], b_rhs = rhs[1];
  double r2, lap_a, lap_b;
  foreach(){
    r2 = sq(x) + sq(y); lap_a = 0; lap_b = 0;
    foreach_dimension(){
      lap_a += (a[1] - 2.0*a[] + a[-1])/sq(Delta);
      lap_b += (b[1] - 2.0*b[] + b[-1])/sq(Delta);
    }
    a_rhs[] = 0.5*ka*r2*b[] - (0.5/Mass)*lap_b;
    b_rhs[] =-0.5*ka*r2*a[] + (0.5/Mass)*lap_a;
    // a_rhs[] = lap_a;
    // b_rhs[] = lap_b;
  }
}

a[right] = dirichlet(0);
b[right] = dirichlet(0);
modulus[right] = dirichlet(0);
a[left] = dirichlet(0);
b[left] = dirichlet(0);
modulus[left] = dirichlet(0);
a[top] = dirichlet(0);
b[top] = dirichlet(0);
modulus[top] = dirichlet(0);
a[bottom] = dirichlet(0);
b[bottom] = dirichlet(0);
modulus[bottom] = dirichlet(0);
/**
## Parameters

We change the size of the domain `L0` and set the tolerance of the
implicit diffusion solver. */

int main()
{

  N=256;
  init_grid (N);
  size (1);
  origin(-L0/2, -L0/2, -L0/2);

  periodic (right);
  periodic (top);
  // fp = popen("ppm2mp4 movie.mp4", "w");

  int order = 2;
  dt = 1.0e-6;
  // the initial condition
  foreach() {
    double r2 = sq(x) + sq(y);
    a[] = A0 * exp(-r2 / (2.0*w0*w0));
    b[] = 0.0;
    modulus[] = sq(a[]) + sq(b[]);
  }
  boundary ({a, b, modulus});
  vtk_file();

  for (t = 0; t <= 2.; t += dt) {
    printf("t=%g", t);
    dt = dtnext(1.0e-6);
    printf(" after dt ");
    runge_kutta ({a, b}, t, dt, du, order);
    printf(" after RK ");
    boundary ({a, b, modulus});
    printf(" after BC ");
    foreach()
      modulus[] = sq(a[]) + sq(b[]);
    printf(" after modulus ");
    vtk_file();
    printf(" after vtk end t=%g\n", t);
  }
}


/**
## Outputs

Here we create an mpeg animation of the $C_1$ concentration. The
`spread` parameter sets the color scale to $\pm$ twice the standard
deviation. */




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



/**
Output*/
// #include "output_fields/output_vtu_foreach.h"
// void vtk_file ()
// {
//   int nf = iteration;
//   scalar l[];
//   scalar modulus[];
//   foreach(){
//     l[] = level;
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
//   iteration++;
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
