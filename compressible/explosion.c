/**
# Two- and three-dimensional explosions

We solve the Euler equations for a compressible gas. We also need to
compute volume fractions for the initial condition. */

#include "compressible.h"
#include "fractions.h"
#include "vtk.h"
#if dimension == 2
# define LEVEL 8
#define MINLEVEL 4
#define MAXLEVEL 10
#else // 3D
# define LEVEL 6
#endif

static int iteration=0;

int main() {

    /**
    We make boundary conditions free outflow. */

    foreach_dimension() {
        w.n[right] = neumann(0);
        w.n[left]  = neumann(0);
////        w.n[top] = dirichlet(0);
////        w.t[top] = dirichlet(0);
////        w.n[bottom]  = dirichlet(0);
////        w.t[bottom]  = dirichlet(0);
//
    }
//    bid circle;
//    w.n[circle] = dirichlet(0);

//    mask (sq(x - 0.5) + sq(y - 0.5) < sq(0.5) ? circle : none);

    /**
    The domain spans $[-1:1]\times[-1:1]\times[-1:1]$. */

    origin (-1, -1, -1);
    size (2.); //size of the problem!
    init_grid (1 << LEVEL);
    run();
}

/**
Initial conditions come from Toro's book (Riemann Solvers and
Numerical Methods for Fluid Dynamics, 3rd Edition, Springer Ed.)
Chapter 17 section 17.1.1 are given in terms of density ($\rho$),
pression ($p$), velocity ($u$) both at the left and right side of the
discontinuity placed at $R=0.4$. */
struct MyCoord
{
    double x;
    double y;
    double z;
//    struct MyCoord (double x_, double y_, double z_){
//        x = x_;
//        y = y_;
//        z = z_;
//    }
};
scalar f[], f1[], f2[], f3[], f4[], f5[], f6[], f7[], f8[], f9[];

void bubble(const struct MyCoord c0, const double R, scalar ftmp) {
    double x0=c0.x, y0=c0.y, z0=c0.z;
    printf ("t = %g c0=%g,%g,%g\n", t, c0.x, c0.y, c0.z);
    fraction (ftmp, sq(x-x0) + sq(y-y0) + sq(z-z0) - sq(R));
    printf ("after\n");
}

event init (t = 0)
{
    double R = 0.1 ;
    double rhoL = 1., rhoR = 0.125 ;
    double pL = 1.0,  pR = 0.1 ;
//    struct MyCoord c0 = {0,0,0};
    /**
    To define an accurate (i.e. symmetrical) initial sphere of rayon
    $R$, we compute the volume fraction corresponding to the spherical
    interface. */


//    bubble(c0, R);
//    f1[] = ftmp[];
//    c0.x= 0.5;
//    bubble(c0, R);
//    f2[] = ftmp[];
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
//        f[] = f1[];

        f[] = min(f1[],min(f2[],min(f3[],min(f4[],min(f5[],min(f6[],min(f7[],min(f8[],f9[]))))))));
//        printf ("i = %d t = %g f1 = %g f2 = %g f = %g\n", i, t, f1[], f2[], f[]);
    }
    /**
    Left and right initial states for $\rho$, $\mathbf{w}$ and energy
    $E = \rho \mathbf{u}^2/2 + p/(\gamma-1)$. */

    foreach() {
        rho[] = rhoR*f[] + rhoL*(1. - f[]);
        foreach_dimension()
        w.x[] = 0.;
        E[] = (pR*f[] + pL*(1. - f[]))/(gammao - 1.);
    }

    theta = 1.3; // tune limiting from the default minmod
}
//
//#include "view.h"
//event image (t = end) {
//    clear();
//    draw_vof ("f");
//    box();
//    save ("image.ppm");
//}

event images (t+= 4./300.) {
//    output_ppm (rho, linear = true);

    scalar l[];
    foreach()
    l[] = level;
    static FILE * fp = fopen ("grid.ppm", "w");
    output_ppm (l, fp, min = 0, max = LEVEL);

//    static FILE * fprho = fopen ("out", "w");
//    output_ppm (rho, fprho, min = 0, max = 1);
//    char name[80];
//    printf("iter=%d\n", iteration);
//    sprintf(name, "list.vtk.%d", iteration++);
//    printf("name %s", name);
//    FILE * fplist = fopen (name, "w");
//    output_vtk ({rho, w}, N, fplist,  true);

}


event movie (t += 0.01; t <= 1.3) {
    static FILE * fprho = popen ("ppm2mpeg > rho.mpg", "w");
    output_ppm (rho, fprho, linear = true);
    scalar wn[];
    foreach(){
        wn[]=0;
        foreach_dimension(){
            wn[]+=w.x[]*w.x[];
        }
        wn[]=sqrt(wn[]);
    }
    static FILE * fpw = popen ("ppm2mpeg > normvel.mpg", "w");
    output_ppm (wn, fpw, linear = true);
}
event end (t = 1.3) {
    //printf ("i = %d t = %g\n", i, t);
}
/**
On trees, we adapt the mesh by controlling the error on the density
field. */

#if TREE
event adapt (i++) {
  adapt_wavelet ({rho}, (double[]){5e-3}, LEVEL + 1);
}
#endif

/**
## Results

Results are presented in terms of $\rho$ and normal velocity $u_n$ for
Cartesian (7 levels in 2D and 6 levels in 3D) and adaptive (8 levels
in 2D and 7 levels in 3D) computations. The numerical results compare
very well with Toro's numerical experiments.

~~~gnuplot Radial density profile
set xrange [0:1]
set xlabel 'r'

set term PNG enhanced font ",10"
set output 'rho.png'
set ylabel 'rho'
plot './cout' u 1:2 w p pt 7 ps 0.2 t '2D Cartesian', \
     './out' u 1:2 w p pt 7 ps 0.2 t '2D Adaptive', \
     '../explosion3D/out' u 1:2 w p pt 7 ps 0.2 t '3D Cartesian', \
     '../explosion.3D/out' u 1:2 w p pt 7 ps 0.2 t '3D Adaptive'
~~~

~~~gnuplot Normal velocity
set output 'velocity.png'
set ylabel 'Normal velocity'
plot './cout' u 1:3 w p pt 7 ps 0.2 t '2D Cartesian',		  \
     './out' u 1:3 w p pt 7 ps 0.2 t '2D Adaptive',		  \
     '../explosion3D/out' u 1:3 w p pt 7 ps 0.2 t '3D Cartesian', \
     '../explosion.3D/out' u 1:3 w p pt 7 ps 0.2 t '3D Adaptive'
~~~

*/

//event print (t = 0.25)
//{
//
//    /**
//    At $t=0.25$ we output the values of $\rho$ and the normal velocity
//    $\mathbf{u}_n$ as functions of the radial coordinate. */
//
//    foreach() {
//        double r = sqrt(sq(x) + sq(y) + sq(z));
//        double wn = (w.x[]*x + w.y[]*y + w.z[]*z)/r;
//        printf ("%g %g %g\n", r, rho[], wn/rho[]);
//    }
//
//    /**
//    For reference we also output a cross-section at $y=0$. */
//
//    for (double x = 0; x <= 1; x += 1e-2)
//        fprintf (stderr, "%g %.4f %.4f\n", x,
//                 interpolate (rho, x, 0.),
//                 interpolate (w.x, x, 0.));
//}