//We want to reproduce the geometry from this page: https://en.wikipedia.org/wiki/Constructive_solid_geometry
//
//Utilisation of boolean operation for geometric construction
//The goal is to define a geometry by using several boolean operation. This geometry is describe here

#include "grid/octree.h"
#include "utils.h"
#include "fractions.h"
#include "view.h"
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

double cubeF(double x, double y, double z, coord center, double size) {
//    Our cube is defined as the intersection of 6 orthogonal planes. We define first 2 planes,  and .
//
//    Then we define P1 has:

    double P1_Plus = x - size/2. + center.x;
    double P1_Minus = x + size/2. + center.x;

    double P1 = max (P1_Plus, -P1_Minus);
//    We apply the same process to obtain P2 and P3

    double P2_Plus = y - size/2. + center.y;
    double P2_Minus = y + size/2. + center.y;

    double P2 = max (P2_Plus, -P2_Minus);

    double P3_Plus = z - size/2. + center.z;
    double P3_Minus = z + size/2. + center.z;

    double P3 = max (P3_Plus, -P3_Minus);
//    At the end, our cube is:


    double c = max ( P1, max (P2, P3) );

    return c;
}
//The sphere function will return a sphere, centered on “center”, of radius “radius”

double sphere(double x, double y, double z, coord center, double radius) {
    return ( sq(x - center.x) + sq (y - center.y) + sq (z - center.z)
             - sq (radius));
}
//We will generate the CSG geometry examples. This geometry is define by using a cube (C); a sphere (S) and 3 cylinders oriented allong ,  and (respectively ,  and ).


double geometry(double x, double y, double z) {

    coord center;
    foreach_dimension()
    center.x = 0;
//    We define the sphere and the cube

    double s = sphere (x, y, z, center, 0.25);
    double c = cubeF (x, y, z, center, 0.38);
//    sIc is the intersection between the sphere and the cube.

    double sIc = max (s, c);
//    We define the 3 cylinders along the ,  and  axis.

    double cylinder1 = sq(x) + sq(y) - sq(0.12);
    double cylinder2 = sq(z) + sq(y) - sq(0.12);
    double cylinder3 = sq(x) + sq(z) - sq(0.12);
//    We use an intermediate object for the union (to decompose the process)

    double cylinderInter = min(cylinder1, cylinder2);
    double cylinderUnion = min(cylinderInter, cylinder3);
//    Finally, we define the geometry we will return.

    double geom = max(sIc, -cylinderUnion);

    return geom;
}
//We have define all the important function for the geometry generation. We can now implemet that into a basilisk mesh.

int main() {
//    We shift the origin so that the simulation domain will be centered on

    origin(-0.5, -0.5, -0.5);
//    We initialise the grid with 7 levels.

    init_grid(1<<7);
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

            view (fov = 13.0359, quat = {-0.353553,0.353553,0.146447,0.853553},
                  tx = 0, ty = 0, bg = {0.3,0.4,0.6}, width = 800, height = 800, samples = 4);
    draw_vof("f", edges = true);
    save ("vof.png");
}