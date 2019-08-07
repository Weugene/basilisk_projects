//An example of 2D viscous flow around a simple solid boundary. Fluid is injected to the left of a channel bounded by solid walls with a slip boundary condition. The Reynolds number is set to 160.

#include "navier-stokes/centered.h"
//#include "grid/quadtree.h"
#define LEVEL 9
#define MINLEVEL 4
#define MAXLEVEL 12

//The domain is eight units long, centered vertically.

int main() {
    L0 = 8.;
    origin (-0.5, -L0/2.);
    //N = 512;
    init_grid(1<<LEVEL);
    run();
}
//The fluid is injected on the left boundary with a unit velocity. The tracer is injected in the lower-half of the left boundary. An outflow condition is used on the right boundary.

u.n[left]  = dirichlet(1.);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);


u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);
//We add a new boundary condition for the cylinder. The tangential velocity on the cylinder is set to zero.

bid cylinder;
u.t[cylinder] = dirichlet(0.);

event init (t = 0) {
//    To make a long channel, we set the top boundary for  and the bottom boundary for . The cylinder has a radius of 0.0625.
    printf("begin init");
    mask (y >  0.5 ? top :
          y < -0.5 ? bottom :
          sq(x) + sq(y) < sq(0.0625) ? cylinder :
          sq(x - 0.2*L0) + sq(y) < sq(0.0625) ? cylinder :
          sq(x - 0.4*L0) + sq(y) < sq(0.0625) ? cylinder :
          sq(x - 0.6*L0) + sq(y) < sq(0.0625) ? cylinder :
          sq(x - 0.8*L0) + sq(y) < sq(0.0625) ? cylinder :
          none);
//    We set a constant viscosity corresponding to a Reynolds number of 160, based on the cylinder diameter (0.125) and the inflow velocity (1). We also set the initial velocity field and tracer concentration.

    const face vector muc[] = {0.00078125,0.00078125};
    mu = muc;
    foreach() {
        u.x[] = 1.;
    }
    printf("end init");
}
//We check the number of iterations of the Poisson and viscous problems.

event logfile (i++)
    fprintf (stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);



event movies (t += 0.1; t <= 15.) {
    static FILE * fp = popen ("ppm2mpeg > vort.mpg", "w");
    scalar vorticity[];
    foreach()
    vorticity[] = (u.x[0,1] - u.x[0,-1] - u.y[1,0] + u.y[-1,0])/(2.*Delta);
    boundary ({vorticity});

    output_ppm (vorticity, fp, box = {{-0.5,-0.5},{7.5,0.5}},
        min = -10, max = 10, linear = true);
}

event images (t += 0.1) {
//    output_ppm (rho, linear = true);

    scalar l[];
    foreach()
    l[] = level;
    static FILE * fp = fopen ("grid.ppm", "w");
    output_ppm (l, fp, min = 0, max = MAXLEVEL);

    static FILE * fpv = fopen ("vort.ppm", "w");
    scalar vorticity[];
    foreach()
        vorticity[] = (u.x[0,1] - u.x[0,-1] - u.y[1,0] + u.y[-1,0])/(2.*Delta);
    boundary ({vorticity});
    output_ppm (vorticity, fpv, linear=true);

    static FILE * fpp = fopen ("p.ppm", "w");
    output_ppm (p, fpp, linear=true);
//    char name[80];
//    printf("iter=%d\n", iteration);
//    sprintf(name, "list.vtk.%d", iteration++);
//    printf("name %s", name);
//    FILE * fplist = fopen (name, "w");
//    output_vtk ({p, u.x, u.y}, N, fplist,  true);

}

#if TREE
event adapt (i++) {
  adapt_wavelet ({p,u}, (double[]){1e-2,1e-3,1e-3}, MAXLEVEL);
}
#endif