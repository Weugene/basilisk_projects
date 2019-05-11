/**
# Two-phase flow around RV fiber

This is an improved version of the [famous Gerris
example](http://gerris.dalembert.upmc.fr/gerris/examples/examples/fiber.html),
illustrating the combination of complex solid boundaries, air-water
turbulent flows and reduced gravity approach.

We use the centered Navier--Stokes solver, two-phase flow and the
momentum-conserving option. Note that the momentum-conserving option
is crucial to obtain stable solutions for this air-water density ratio
configuration. */

//#include "grid/octree.h"
#include "grid/quadtree.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"

#include "view.h"
#include "output_vtu_foreach.h"
//#include "vtk.h"
//#include "output_vtk.h"


int MAXLEVEL = 10;
int LEVEL = 6;
float file_step = 1000;
float file0 = 0;
char file_template[80] = "dump-%g";
# define MAXITER 1000000


/**
## Main function
 */

scalar fiber[], f0[];

int main (int argc, char * argv[])
{

    L0 =4;
    printf("write file template: dump-%%g, first file num, step...\n");
    if (argc > 1) {
        snprintf(file_template, 80, "%s", argv[1]);
    }
    if (argc > 2) {
        file0 = atof(argv[2]);
    }
    if (argc > 3) {
        file_step = atof(argv[3]);
    }
    printf("template: %s\n", file_template);
    printf("first file num = %g\n", file0);
    printf("step = %g\n", file_step);

    init_grid (64);

    size (L0);
    origin (-L0/2.,-L0/2.,-L0/2.);
    run();
}

/**
## Initial conditions

We can optionally restart, otherwise we open the STL file and
initialize the corresponding fraction. We also initialize the `f0`
field used for the inflow condition and set the initial water level
and velocity field. */

event init (i = 0) {
    int k=0;
    while(1) {
        char name[80];
        sprintf(name, file_template, file0 + file_step*k);
        printf("dump name=%s\n", name);
        if (restore(file = name)) {
//            printf("begin of reading...\n");
//            char namevtk[80];
//            sprintf(namevtk, "list.vtk.%d", k);
//            printf("namevtk = %s\n", namevtk);
//            FILE *fplist = fopen(namevtk, "w");
            scalar l[];
            foreach()
                l[] = level;
//            output_vtk({f, fiber, u.x, u.y, rhov, p, pf, l}, 1024, fplist, false);
//            output_pvtu_bin (scalar * list, vector * vlist, int n, FILE * fp, char * subname)
//            printf("%s generated\n", name);
//
            char name[80], subname[80];
            FILE *fp;
            sprintf(name, "out_%4.4d_n%3.3d.vtu", k, 0);
            fp = fopen(name, "w");
            output_vtu_bin_foreach((scalar *) {f, l, fiber, p, rhov}, (vector *) {u}, 64, fp, false);
            fclose(fp);

            sprintf(name, "out_%4.4d.pvtu", k);
            sprintf(subname, "out_%4.4d", k);
            fp = fopen(name, "w");
            output_pvtu_bin((scalar *) {f, l, fiber, p, rhov}, (vector *) {u}, 64, fp, subname);
            fclose(fp);

        }else{
            printf("end of reading\n");
            return 888;
            break;
        }
        k++;
    }
    printf("init is read");
}



