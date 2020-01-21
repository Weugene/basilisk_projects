/**
We use the centered Navier--Stokes solver, two-phase flow and the
momentum-conserving option. Note that the momentum-conserving option
is crucial to obtain stable solutions for this air-water density ratio
configuration. */

#define BRINKMAN_PENALIZATION 1
#define DEBUG_BRINKMAN_PENALIZATION 1
#define DEBUG_MINMAXVALUES 1

#undef SEPS
#define SEPS 1e-30

//#include "grid/octree.h"
#include "grid/quadtree.h"
#include "../src_local/centered-weugene.h"// u, uf, p,  pf, g, mu=0, alphav=1, a=0(acceleration), rho=1,
// stokes=false (No omit div(u by u)
// incompressible NS.  generic time loop, a CFL-limited timestep,
// the Bell-Collela-Glaz (BCG) advection scheme ( 2nd-order, unsplit, upwind scheme)
// and the implicit viscosity solver
/**
We have two phases e.g. air and water. For large viscosity and density
ratios, the harmonic mean for the viscosity tends to work better than
the default arithmetic mean. We "overload" the default by defining the
*mu()* macro before including the code for two phases. */
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h" //advection of tracer f. defined rhov, alphav. if FILTERED is defined, then sf - filtered values of f.
#include "navier-stokes/conserving.h" //???should add?
#include "../src_local/output_vtu_foreach.h"

#include "tension.h"
#include "view.h"
//#include "distance.h"
#include "reduced.h"

/**
On supercomputers we need to control the maximum runtime and we check
performances. */

//#include "maxruntime.h"
//#include "navier-stokes/perfs.h"

/**
## Importing the geometry

This function computes the solid fraction given a pointer to an STL
file, a tolerance (maximum relative error on distance) and a
maximum level. */

//void fraction_from_stl (scalar f, FILE * fp, double eps, int maxlevel)
//{
//
//    /**
//    We read the STL file and compute the bounding box of the model. */
//
//    coord * p = input_stl (fp);
//    coord min, max;
//    bounding_box (p, &min, &max);
//    double maxl = -HUGE;
//    foreach_dimension()
//    if (max.x - min.x > maxl)
//        maxl = max.x - min.x;
//
//    /**
//    We initialize the distance field on the coarse initial mesh and
//    refine it adaptively until the threshold error (on distance) is
//    reached. */
//
//    scalar d[];
//    distance (d, p);
//    while (adapt_wavelet ({d}, (double[]){eps*maxl}, maxlevel, 5).nf);
//
//    /**
//    We also compute the volume fraction from the distance field. We
//    first construct a vertex field interpolated from the centered field
//    and then call the appropriate VOF functions. */
//
//    vertex scalar phi[];
//    foreach_vertex()
//    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1] +
//             d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.; //left bottom vertex
//    fractions (phi, f);
//}


int maxlevel = 10;
int minlevel = 5;
int LEVEL = 10;
double limMin, limMax;
/**
The density ratio is 1000 and the dynamic viscosity ratio 100. */
#define RHOR 1000.
#define MUR 100.
#define feps 0.001
#define uemax 0.01
/**
* We consider two bubbles studied by Cano-Lozano et al, 2016. */
# define RE 100.25
# define CA 0.1


#define RandMinMax(limMin, limMax) (0.5*((limMin)+(limMax)) + 0.5*noise()*((limMax) - (limMin)))
#if dimension>2
# define mynorm(v) (sqrt(sq(v.x) + sq(v.y) + sq(v.z)))
#else
# define mynorm(v) (sqrt(sq(v.x) + sq(v.y)))
#endif
/**
 If the distance function is positive, then the value of $f=1$, otherwise it becomes 0.
 For example, at some point $f=1$ means that a point consist of pure liquid, and
 $f=0$ means that there is only gas in a point, no liquid.
 */

double cubeF(double x, double y, double z, coord center, double size) {
//    Our cube is defined as the intersection of 6 orthogonal planes. We define first 2 planes,  and .
//
//    Then we define P1 has:
    double P1_Plus = x - size/2. + center.x;
    double P1_Minus = x + size/2. + center.x;
    double P1 = max (P1_Plus, - P1_Minus);
//    We apply the same process to obtain P2 and P3
    double P2_Plus = y - size/2. + center.y;
    double P2_Minus = y + size/2. + center.y;
    double P2 = max (P2_Plus, -P2_Minus);
#if dimension>2
    double P3_Plus = z - size/2. + center.z;
    double P3_Minus = z + size/2. + center.z;
    double P3 = max (P3_Plus, -P3_Minus);
    double c = max ( P1, max (P2, P3) );
#else
    double c = max ( P1, P2 );
#endif
    return -c;
}
double sphere(double x, double y, double z, coord center, double radius) {
    return ( sq(x - center.x) + sq (y - center.y) + sq (z - center.z)
             - sq (radius));
}

double bubbles (double x, double y, double z)
{
    const int ns = 10;
    coord mypoint = {x,y,z};
    coord pnt_dist;
    coord centers[ns];
    double R[ns], dist=1.1*sqrt(sq(limMax - limMin)+sq(L0)
#if dimension>2
            +sq(L0)
#endif
    )/(ns+2.0);;

    srand (0);
    int iter = 0, i = 0;
    while(i<ns){
        R[i] = 0.01*L0*(fabs(noise())+1.0);
        centers[i].x = RandMinMax( limMin + R[i], limMax - R[i]);
        centers[i].y = RandMinMax(-L0/2.0, L0/2.0);
#if dimension>2
        centers[i].z = RandMinMax(-0.5*L0 + R[i], 0.5*L0 - R[i]);
#endif
        for (int j = 1; j < i; j++) {
            foreach_dimension()
            pnt_dist.x = centers[i].x - centers[j].x;
            if (mynorm(pnt_dist) < dist) {
                i--;
                break;
            };
        }
        i++;
        iter++;
        if (iter>100*ns) exit(137);
    }

    for (int i = 0; i < ns; i++){
        R[i] = 0.02*L0*(fabs(noise())+1.0);
        centers[i].x = RandMinMax( limMin + R[i], limMax - R[i]);
        centers[i].y = RandMinMax(-0.5*L0, 0.5*L0);
//        centers[i].y = RandMinMax(-7.0*R[i], 7.0*R[i]);
#if dimension>2
        centers[i].z = RandMinMax(-0.5*L0 + R[i], 0.5*L0 - R[i]);
#endif
        for (int j = 0; j < i; j++) {
            foreach_dimension() pnt_dist.x = centers[i].x - centers[j].x;
            if (mynorm(pnt_dist) < dist) {
                i--;
                break;
            };
        }
    }

    double phi = HUGE;
    for (int i = 0; i < ns; i++) {
        foreach_dimension() pnt_dist.x = mypoint.x - centers[i].x;
        phi = min(phi, (mynorm(pnt_dist) - R[i]));
    }
//    return min(phi, -L0/8.0 - x);// with front
    return phi; // no front
}

coord center={0,0,0};
const int Ncyl=7;
double size_box;
double R;//0.5
double dist;

double cylindersOy(double x, double y, double z, double Ox, int Ncyl) {
    int icyl;
//    We define the 2*Ncyl cylinders along the ,  and  axis.
    double cylinderX[Ncyl], tmp;
#if dimension>2
    double  cylinderZ[Ncyl];
#endif
    for (icyl=0; icyl<Ncyl; icyl++){
//        tmp = dist*(icyl-(Ncyl-1.0)/2.0);
        tmp = dist*icyl + center.y - 0.5*size_box + R;
//        cylinderX[icyl] = sq(R)- sq(y+2.0*R) - sq(x-tmp);
        cylinderX[icyl] = sq(R)- sq(x + 2.0*R - Ox) - sq(y - tmp);
#if dimension>2
        cylinderZ[icyl] = sq(R) - sq(y) - sq(z-tmp) ;
#endif
    }
    //    We use an intermediate object for the union (to decompose the process)
    double geom = cylinderX[0];
    for (icyl=1; icyl<Ncyl; icyl++){
        geom = max(cylinderX[icyl], geom);
    }
#if dimension>2
    for (icyl=0; icyl<Ncyl; icyl++){
        geom = max(cylinderZ[icyl], geom);
    }
#endif
    return geom;
}

double geometry(double x, double y, double z) {
    const int Nlayer=4;
    int icyl;
    double dist = 0.5*L0/Nlayer;
    double geom = cylindersOy(x, y, z, 0.0, Ncyl);

    for (icyl=1; icyl<Nlayer; icyl++){
        geom = max(cylindersOy(x, y+dist/2.*(icyl%2), z, icyl*dist, Ncyl + icyl%2), geom);
    }
    double c = cubeF(x, y, z, center, size_box);
    return min(geom,c);
}
/**
## Main function
We need additional (fraction) fields for the ship geometry and for the
(inflow) boundary condition. */

scalar fs[];

int main (int argc, char * argv[])
{
    maxruntime (&argc, argv);
    if (argc > 1) {
        LEVEL = atoi(argv[1]); //convert from string to int
    }
    periodic (top);
#if dimension>2
    periodic (back);
#endif
    init_grid (64);
    rho1 = 1.;// water
    rho2 = 1./RHOR; // air
    mu1 = 1./RE;
    mu2 = 1./(MUR*RE);
    f.sigma = 1./CA;
    eta_s = 1e-15;
    /**
    The length of the ship is unity and the domain is five times
    larger. We change the origin so that the ship is not too close to
    the inflow. */
    stokes = true;
    size (4);
    origin (-L0/2.,-L0/2.,-L0/2.);
    limMin = -0.5*L0;
    limMax = -0.25*L0;
    /**
  We reduce the tolerance on the divergence of the flow. This is
  important to minimise mass conservation errors for these simulations
  which are very long. */

    TOLERANCE = 1e-8;
    for (scalar s in {f,fs}) {
        s.refine = s.prolongation = fraction_refine;
    }

    run();
}

/**
## Boundary conditions*/

u.n[left] = neumann(0);
p[left]   = dirichlet(1);
pf[left]  = dirichlet(1);
f[left]   = dirichlet(1); //liquid

u.n[right] = neumann(0);
p[right]   = dirichlet(0);
pf[right]  = dirichlet(0);
f[right]   = neumann(0);//liquid

/**
## Initial conditions

We can optionally restart, otherwise we open the STL file and
initialize the corresponding fraction. We also initialize the `f`
field used for the inflow condition and set the initial water level
and velocity field. */

event init (t = 0) {
    size_box = 0.8*L0;
    R = size_box*0.2/(max(Ncyl-1,1));
    dist = (size_box-2.0*R)/(max(Ncyl-1,1));

    if (!restore (file = "restart")) {
//        from stl file
//        FILE * fp = fopen ("fiber.stl", "r");
//        fraction_from_stl (fs, fp, 5e-4, maxlevel);
//        fclose (fp);
        int iteration = 0;
        do {
            iteration++;
            fraction(fs, geometry (x, y, z));
            fraction(f, bubbles (x, y, z));
        }while (adapt_wavelet({fs, f}, (double []){feps, feps},
                              maxlevel = maxlevel, 2).nf != 0 && iteration <= 10);
        if (pi) fprintf(ferr, "initializing of fs f have been finished \n");
        foreach() {
            foreach_dimension() u.x[] = 0.;
        }
        boundary ({f, u.x, u.y});
    }
}

//event movie (t += 0.01; t <= 10)
//{
//#if dimension == 2
//view (tx = 0, width = 2048, height = 2048, samples = 4, bg = {0.3,0.4,0.6});
//    clear();
//    draw_vof ("f", edges = true);
//    draw_vof("fiber", edges = true);
//    box ();
//    save ("vof2D.mp4");
//#else // 3D
//clear();
//view (fov = 40,
//        quat = {0.5,0.1,0.2,0.8},
//        tx = 0, ty = 0,
//        width = 1600, height = 1600);
//draw_vof ("fiber");
//draw_vof ("f", linear = true);
//save ("movie.mp4");
//#endif
//
//}

//Output
event vtk_file (t += 0.01; t<100)
{
    scalar l[]; foreach() l[] = level;
    char subname[80]; sprintf(subname, "2phase");
    output_vtu_MPI( (scalar *) {l, f, fiber, p}, (vector *) {u}, subname,  L0/pow(2, minlevel));
}

#if DUMP
    event snapshot (t += 0.01) {
      char name[80];
      sprintf (name, "dump-%d", i);
      fprintf(ferr, "dump: %s, t=%g, dt=%g, i=%d", name, t, dt, i);
      dump (file = name);
    }
#endif

/**
## Mesh adaptation

This computation is only feasible thanks to mesh adaptation, based
both on volume fraction and velocity accuracy. */
#define ADAPT_SCALARS {f, fs, u}
#define ADAPT_EPS_SCALARS {feps, feps, uemax, uemax, uemax}
event adapt (i++){
    double eps_arr[] = ADAPT_EPS_SCALARS;
    //    MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = 5);
    //    calc_solid(fs, n_sol, target_U);
}




