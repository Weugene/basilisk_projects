//We want to reproduce the geometry from this page: https://en.wikipedia.org/wiki/Constructive_solid_geometry
//
//Utilisation of boolean operation for geometric construction
//The goal is to define a geometry by using several boolean operation. This geometry is describe here

#include "grid/octree.h"
#include "utils.h"
#include "fractions.h"
#include "view.h"

#include "navier-stokes/centered.h"

#define LEVEL 9
#define MINLEVEL 4
#define MAXLEVEL 12
//Definition of the boolean operation
//Between 2 volumes (or surface in 2D), we can define a union, an intersection and a substraction of the volume.
//
//Let use A and B, 2 volumes of the space, such that  and . Then:
//
//
//
//Pay attention, in Basilisk, the surface expressions have to be continuous function.
//
//For example, a sphere is define with the following function:
//
//
//In the sphere,  is negative, allong the interface,  is equal to zero and outside of the sphere,  is positive. The function is continuous.
//
//In basilisk, once you gave  to the macros “fraction”, the code will generate an interface when  is equal to zero.
//
//Definition of the geometric interface
//        We define first 2 basic geometric functions, a shpere and a cube.
//
//The cube is centered on “center”, and has a lenght of “size”.


//We will generate the CSG geometry examples. This geometry is define by using a cube (C); a sphere (S) and 3 cylinders oriented allong ,  and (respectively ,  and ).


double geometry(double x, double y, double z) {
    const int Ncyl=10;
    double R = min(0.05,0.4/Ncyl);//0.5
    double dist = (1.0-2.0*R)/(Ncyl-1.0);

    int icyl;
//    We define the 2*Ncyl cylinders along the ,  and  axis.
    double cylinderX[Ncyl], cylinderY[Ncyl], tmp;
    for (icyl=0; icyl<Ncyl; icyl++){
//        tmp = dist*(icyl-(Ncyl-1.0)/2.0);
        tmp = dist*icyl + R - 0.5;
        cylinderX[icyl] = sq(z+2.0*R+dist) + sq(y-tmp) - sq(R);
        cylinderY[icyl] = sq(z) + sq(x-tmp) - sq(R);
    }
    //    We use an intermediate object for the union (to decompose the process)
    double geom = min(cylinderX[0], cylinderY[0]);
    for (icyl=1; icyl<Ncyl; icyl++){
        geom = min(min(cylinderX[icyl], cylinderY[icyl]),geom);
    }

    return geom;
}
//We have define all the important function for the geometry generation. We can now implement that into a basilisk mesh.

int main() {
//    We shift the origin so that the simulation domain will be centered on

    origin(-0.5, -0.5, -0.5);
//    We initialise the grid with 7 levels.

    init_grid(1<<LEVEL);
//    We create the scalar field that will get the geometry.

    scalar f[];
//    To have a refine mesh on the geometry, we initialise an iterative compt

    int iteration = 0;
    do {
        iteration++;
        fraction(f, geometry (x, y, z));
    }while (adapt_wavelet({f}, (double []){0.2},
                          maxlevel = 9, 2).nf != 0 && iteration <= 10);
//    The geometry obtain with basilisk can be observe by using bview. The generated picture is:

//    Reconstructed VOF surface.
//            Reconstructed VOF surface.

    view (fov = 50, quat = {-0.4,0.4,0.1,0.8}, tx = 0, ty = 0, bg = {0.3,0.4,0.6}, width = 800, height = 800, samples = 4);
    draw_vof("f", edges = true);
    save ("vof.png");
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
    mask (geometry (x, y, z) < 0 ?cylinder:
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
#if dimension == 2
    scalar omega[];
    vorticity (u, omega);
      view (tx = -0.5);
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