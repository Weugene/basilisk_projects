/**
# Stokes flow through a complex 3D porous medium

This is the 3D equivalent of the [2D porous medium test
case](/src/test/porous.c).

The medium is periodic and described using embedded boundaries.

This tests mainly the robustness of the representation of embedded
boundaries and the convergence of the viscous and Poisson
solvers. */
#define RELATIVE_RESIDUAL
#define FILTERED
#define DEBUG_MINMAXVALUES
#define DEBUG_OUTPUT_VTU_MPI
#define EPS_MAXA 2
#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
#include "view.h"
#include "navier-stokes/perfs.h"
#include "../src_local/output_vtu_foreach.h"
#include "../src_local/utils-weugene.h"
/**
We also need to compute distance functions (to describe the solid geometry
), use visualisation functions. */

#include "distance.h"
/**
We will vary the maximum level of refinement, starting from 5. */

int minlevel = 4;
int maxlevel = 10;
double cseps = 1e-3, ueps = 1e-3;
double RE = 42, FROUDE = 0.7;
const face vector muv_[] = {0.1,0.1,0.1};
face vector muv[];
/**
The porous medium is defined by the union of a random collection of
spheres. The number of spheres `ns` can be varied to vary the porosity. */

void porous (scalar cs, face vector fs)
{
    int ns = 700;
    coord pc[ns];
    double R[ns];
    srand (0);
    for (int i = 0; i < ns; i++) {
        foreach_dimension()
        pc[i].x = 0.5*noise();
        R[i] = 0.04 + 0.08*fabs(noise());
    }

    /**
    Once we have defined the random centers and radii, we can compute
    the levelset function $\phi$ representing the embedded boundary. */

    vertex scalar phi[];
    foreach_vertex() {
        phi[] = HUGE;

        /**
        Since the medium is periodic, we need to take into account all the
        disk images using periodic symmetries. Note that this means that
        each function evaluation requires 27 times `ns` evaluations of the
        function for a sphere. This is expensive but could be improved a
        lot using a more clever algorithm. */

        for (double xp = -L0; xp <= L0; xp += L0)
            for (double yp = -L0; yp <= L0; yp += L0)
                for (double zp = -L0; zp <= L0; zp += L0)
                    for (int i = 0; i < ns; i++)
                        phi[] = intersection (phi[], (sq(x + xp - pc[i].x) +
                                                      sq(y + yp - pc[i].y) +
                                                      sq(z + zp - pc[i].z) -
                                                      sq(R[i])));
    }
    boundary ({phi});

    fractions (phi, cs, fs);

    /**
    This is necessary to remove degenerate fractions which could cause
    convergence problems. */

    fractions_cleanup (cs, fs);
}

void bone (scalar cs, face vector fs)
{
    double eps = 1e-3;
    FILE * fp = fopen ("cube.stl", "r");
    fprintf (ferr, "STL file is opened... \n");
    /**
   We read the STL file and compute the bounding box of the model. */
    coord * p = input_stl (fp);
    coord min, max;
    bounding_box (p, &min, &max);
    fprintf (ferr, "min = (%g %g %g), max = (%g %g %g) \n", min.x, min.y, min.z, max.x, max.y, max.z);
    double maxl = -HUGE;
    foreach_dimension() if (max.x - min.x > maxl)
        maxl = max.x - min.x;

    /**
    We initialize the distance field on the coarse initial mesh and
    refine it adaptively until the threshold error (on distance) is
    reached. */

    scalar d[];
    int it = 0;
    distance (d, p);
    while (adapt_wavelet ({d}, (double[]){eps*maxl}, maxlevel=maxlevel, minlevel=minlevel).nf && ++it<=10);
    fprintf (ferr, "after distance...\n");
    int no_cells=0;
    foreach (reduction (+:no_cells)){
        no_cells++;
    }
    int nc_max = (int) pow (2, dimension*maxlevel);
    fprintf (ferr, "it= %d number of cells = %d, max number of cells = %d, compression ratio = %g\n", it, no_cells, nc_max, (double) no_cells/nc_max );

    /**
    We also compute the volume fraction from the distance field. We
    first construct a vertex field interpolated from the centered field
    and then call the appropriate VOF functions. */

    vertex scalar phi[];
    foreach_vertex()
    phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1]
             + d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.; //left bottom vertex
    boundary ({phi});
    fractions (phi, cs, fs);
    fprintf (ferr, "STL file is saved in fs... \n");
    fclose (fp);
}

/**
The boundary condition is zero velocity on the embedded boundary. */

u.n[embed] = dirichlet(0);
u.t[embed] = dirichlet(0);
u.r[embed] = dirichlet(0);

/**
The domain is the periodic unit cube centered on the origin. */
int main()
{
    size(1);
//    origin (-0.5, -0.5, -0.5);
    foreach_dimension()
        periodic (right);

    /**
    We turn off the advection term. The choice of the maximum timestep
    and of the tolerance on the Poisson and viscous solves is not
    trivial. This was adjusted by trial and error to minimize (possibly)
    splitting errors and optimize convergence speed. */

    stokes = true;
    DT = 2e-5;
    TOLERANCE = HUGE;
//    NITERMIN = 2;
    N = 1 << 6;

    run();
}

scalar un[];

event init (t = 0) {

    /**
    We define the porous embedded geometry. */

    bone (cs, fs);
//    porous (cs, fs);

    /**
    The gravity vector is aligned with the channel and viscosity is
    unity. */

    const face vector g[] = {1./sq(FROUDE),0.,0.};
    a = g;

//    const face vector muv_[] = {1./RE,1./RE,1./RE};
//    (const) face vector muv = muv_;
//    mu = muv;
//    mu = fm;



    /**
    We initialize the reference velocity. */

    foreach()
        un[] = u.x[];
}

#if 0 // used for debugging only
coord maxpos (scalar s)
{
  coord pmax = {0,0,0};
  double numax = 0.;
  foreach_leaf()
    if (fabs(s[]) > numax) {
      numax = fabs(s[]);
      pmax = (coord){x,y,z};
    }
  return pmax;
}

event dumps (i += 10) {
  scalar nu[], du[], crappy[];
  foreach() {
    nu[] = norm(u);
    du[] = u.x[] - un[];
    double val;
    crappy[] = (embed_flux (point, u.x, mu, &val) != 0.);
  }

#if !_MPI
  coord numax = maxpos (nu);
  fprintf (stderr, "numax: %g %g %g\n", numax.x, numax.y, numax.z);
  numax = maxpos (p);
#if 0
  Point point = locate (numax.x, numax.y, numax.z);
  fprintf (stderr, "pmax: %g %g %g %g %g\n", numax.x, numax.y, numax.z,
	   p[], cs[]);
#endif
  numax = maxpos (du);
  {
    Point point = locate (numax.x, numax.y, numax.z);
    fprintf (stderr, "dumax: %g %g %g\n", numax.x, numax.y, numax.z);
    FILE * fp = fopen ("dumax", "w");
    foreach_neighbor(1) {
      fprintf (fp, "fs %g %g %g %g\n", x - Delta/2., y, z, fs.x[]);
      fprintf (fp, "fs %g %g %g %g\n", x, y - Delta/2., z, fs.y[]);
      fprintf (fp, "fs %g %g %g %g\n", x, y, z - Delta/2., fs.z[]);
      fprintf (fp, "fs %g %g %g %g\n", x + Delta/2., y, z, fs.x[1]);
      fprintf (fp, "fs %g %g %g %g\n", x, y + Delta/2., z, fs.y[0,1]);
      fprintf (fp, "fs %g %g %g %g\n", x, y, z + Delta/2., fs.z[0,0,1]);
      fprintf (fp, "%g %g %g %g\n", x, y, z, cs[]);
    }
    fclose (fp);
  }
#endif

  p.nodump = false;
  dump();
}
#endif

event properties (i++) {
    foreach_face()
        muv.x[] = fm.x[]/RE;
}
/**
We check for a stationary solution. */

event logfile (i++)
{
    double avg = normf(u.x).avg, du = change (u.x, un)/(avg + SEPS);
    fprintf (ferr, "%d %d %d %d %d %d %d %d %.3g %.3g %.3g %.3g %.3g\n",
    maxlevel, i,
    mgp.i, mgp.nrelax, mgp.minlevel,
    mgu.i, mgu.nrelax, mgu.minlevel,
    du, mgp.resa*dt, mgu.resa, statsf(u.x).sum, normf(p).max);

    /**
    If the relative change of the velocity is small enough we stop this
    simulation. */

    if (i > 1 && (avg < 1e-9 || du < 0.01)) {
        /**
        We are interested in the permeability $k$ of the medium, which is
        defined by
        $$
        U = \frac{k}{\mu}\nabla p = \frac{k}{\mu}\rho g
        $$
        with $U$ the average fluid velocity. */

        stats s = statsf (u.x);
        printf ("%d %g\n", maxlevel, s.sum/s.volume);

        /**
        We output fields and dump the simulation. */

        scalar nu[];
        foreach()
                nu[] = norm(u);
                        boundary ({nu});

        view (fov = 32.2073, quat = {-0.309062,0.243301,0.0992085,0.914026},
                tx = 0.0122768, ty = 0.0604286, bg = {1,1,1},
                width = 600, height = 600);
        box();
        draw_vof("cs", "fs", fc = {0.5,0.5,0.5});
        char name[80];
        sprintf (name, "cs-%d.png", maxlevel);
        save (name);

        box();
        isosurface ("u.x", 1e-5, fc = {0.,0.7,0.7});
        sprintf (name, "nu-%d.png", maxlevel);
        save (name);

        view (fov = 19.1765, quat = {0,0,0,1},
                tx = 0.0017678, ty = 0.00157844, bg = {1,1,1},
                width = 600, height = 600);
        squares ("nu", linear = true);
        cells();
        sprintf (name, "cross-%d.png", maxlevel);
        save (name);

        sprintf (name, "dump-%d", maxlevel);
        dump (name);

        /**
        We stop at level 8. */

        if (maxlevel >= 8)
        return 1; /* stop */

        /**
        We refine the converged solution to get the initial guess for the
        finer level. We also reset the embedded fractions to avoid
        interpolation errors on the geometry. */

        maxlevel++;
        #if TREE
        adapt_wavelet ({cs,u}, (double[]){1e-2,4e-6,4e-6,4e-6}, maxlevel);
        #endif
        porous (cs, fs);
        boundary (all); // this is necessary since BCs depend on embedded fractions
        event ("adapt");
    }
}

//event end_timestep (i += 1)
event end_timestep (t += 0.01)
{
    char subname[80]; sprintf (subname, "simple_bone");
    scalar l[];
    foreach() {l[] = level;}
    output_vtu_MPI (subname, t+dt, (scalar *) {cs, p, l}, (vector *) {u} );
}

/**
![Boundary of the porous medium. Fluid is flowing inside this
 volume.](porous3D/cs-8.png)

![Isosurface of the x-component of the velocity.](porous3D/nu-8.png)

![Cross-section at $z = 0$ showing the mesh and norm of the
 velocity.](porous3D/cross-8.png)

~~~gnuplot Permeability as a function of resolution
set xlabel 'Level'
set grid
set ytics format '%.1e'
plot 'out' w lp t ''
~~~

~~~gnuplot Convergence history
set xlabel 'Iterations'
set logscale y
set ytics format '%.0e'
set yrange [1e-8:]
set key center left
plot '../porous3D.ref' u 2:9 w l t '', '' u 2:10 w l t '', \
    '' u 2:11 w l t '', '' u 2:12 w l t '', '' u 2:13 w l t '', \
    'log' u 2:9 w p t 'du', '' u 2:10 w p t 'resp', \
    '' u 2:11 w p t 'resu', '' u 2:12 w p t 'u.x.sum', '' u 2:13 w p t 'p.max'
~~~

## See also

* [Stokes flow past a periodic array of spheres](/src/test/spheres.c)
* [Stokes flow past a periodic array of cylinders](/src/test/cylinders.c)
* [Stokes flow through a complex porous medium](/src/test/porous.c)
*/

/**
## Mesh adaptation

This computation is only feasible thanks to mesh adaptation, based
both on volume fraction and velocity accuracy. */
#define ADAPT_SCALARS {cs, u}
#define ADAPT_EPS_SCALARS {cseps, ueps, ueps, ueps}
event adapt (i++)
{
    double eps_arr[] = ADAPT_EPS_SCALARS;
    MinMaxValues (ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) {u}, (double []){cseps, ueps, ueps, ueps}, maxlevel = maxlevel, minlevel = minlevel);
}

event snapshot (t += 0.01) {
    char name[80];
    sprintf(name, "dump-%04g", t);
    dump (file = name);
}

event stop (t = 17);