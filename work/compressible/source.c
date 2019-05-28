/**
# Two- and three-dimensional explosions

We solve the Euler equations for a compressible gas. We also need to
compute volume fractions for the initial condition. */

#include "compressible.h"
#include "fractions.h"

#if dimension == 2
# define LEVEL 8
#define MINLEVEL 4
#define MAXLEVEL 10
#else // 3D
# define LEVEL 6
#endif

int main() {

    /**
    We make boundary conditions free outflow. */

    foreach_dimension() {
        w.n[right] = neumann(0);
        w.n[left]  = neumann(0);
    }
    origin (-1, -1, -1);
    size (2); //dimension of the problem!
    init_grid (1 << LEVEL);
    run();
}

event init (t = 0)
{
    double R = 0.1 ;
    double rhoL = 1., rhoR = 0.125 ;
    double pL = 1.0,  pR = 0.1 ;
    scalar f[], f1[], f2[], f3[], f4[], f5[], f6[], f7[], f8[], f9[];

    fraction (f1, sq(x) + sq(y) + sq(z) - sq(R));
    fraction (f2, sq(x-0.5) + sq(y) + sq(z) - sq(R));
    fraction (f3, sq(x+0.5) + sq(y) + sq(z) - sq(R));
    fraction (f4, sq(x) + sq(y-0.5) + sq(z) - sq(R));
    fraction (f5, sq(x-0.5) + sq(y-0.5) + sq(z) - sq(R));
    fraction (f6, sq(x+0.5) + sq(y-0.5) + sq(z) - sq(R));
    fraction (f7, sq(x) + sq(y+0.5) + sq(z) - sq(R));
    fraction (f8, sq(x-0.5) + sq(y+0.5) + sq(z) - sq(R));
    fraction (f9, sq(x+0.5) + sq(y+0.5) + sq(z) - sq(R));
    foreach() {
        f[] = min(f1[],min(f2[],min(f3[],min(f4[],min(f5[],min(f6[],min(f7[],min(f8[],f9[]))))))));
    }

    foreach() {
        rho[] = rhoR*f[] + rhoL*(1. - f[]);
        foreach_dimension()
            w.x[] = 0.;
        E[] = (pR*f[] + pL*(1. - f[]))/(gammao - 1.);
    }
    theta = 1.3; // tune limiting from the default minmod
}

event images (t+= 4./300.) {
    output_ppm (rho, linear = true);

    scalar l[];
    foreach()
    l[] = level;
    static FILE * fp = fopen ("grid.ppm", "w");
    output_ppm (l, fp, min = 0, max = LEVEL);
}

event adapt (i++) {
    adapt_wavelet ({rho}, (double []){4e-3}, maxlevel = LEVEL);
}

event end (t = 1.3) {
}
/**
On trees, we adapt the mesh by controlling the error on the density
field. */
#if TREE
event adapt (i++) {
  adapt_wavelet ({rho}, (double[]){5e-3}, LEVEL + 1);
}
#endif
