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
#define EPS_MAXA 2
#include "grid/octree.h"
#include "embed.h"
#include "navier-stokes/centered.h"
//#include "view.h"
#include "navier-stokes/perfs.h"
#include "../src_local/output_vtu_foreach.h"

/**
We also need to compute distance functions (to describe the solid geometry
), use visualisation functions. */

#include "distance.h"

int maxlevel = 8;
int minlevel = 5;
double cseps=1e-3, ueps=1e-3;
#define flux  1e-4 // 4-6 л\мин или 67-100 см^3\сек
#define S 4e-6 //cross section
#define RHO  1050 //плотность 1,050—1,060 г/см³
#define MU 5e-3 //вязкость похоже вот 4-5 мПа×с
#define v_in (flux/S)
void bone (scalar cs, face vector fs){
    double eps = 1e-3;
    FILE * fp = fopen ("cube.stl", "r");
    fprintf (ferr, "STL file is opened... \n");
    /**
   We read the STL file and compute the bounding box of the model. */
    coord * p = input_stl (fp);
    coord min, max;
    bounding_box (p, &min, &max);
    fprintf(ferr, "min = (%g %g %g), max = (%g %g %g) \n", min.x, min.y, min.z, max.x, max.y, max.z);
    double maxl = -HUGE;
    foreach_dimension() if (max.x - min.x > maxl)
        maxl = max.x - min.x;

    /**
    We initialize the distance field on the coarse initial mesh and
    refine it adaptively until the threshold error (on distance) is
    reached. */

    scalar d[];
    int it = 0;
        distance(d, p);
    while (it<10 && adapt_wavelet ({d}, (double[]){eps*maxl}, maxlevel=maxlevel, minlevel=minlevel).nf);
    fprintf(ferr, "after distance...\n");
    int no_cells=0;
    foreach(reduction(+:no_cells)){
        no_cells++;
    }
    int nc_max = (int) pow(2, dimension*maxlevel);
    fprintf(ferr, "it= %d number of cells = %d, max number of cells = %d, compression ratio = %g\n", it, no_cells, nc_max, (double) no_cells/nc_max );

    /**
    We also compute the volume fraction from the distance field. We
    first construct a vertex field interpolated from the centered field
    and then call the appropriate VOF functions. */

    vertex scalar phi[];
    foreach_vertex()
    phi[] = -(d[] + d[-1] + d[0,-1] + d[-1,-1]
             + d[0,0,-1] + d[-1,0,-1] + d[0,-1,-1] + d[-1,-1,-1])/8.; //left bottom vertex
    boundary ({phi});
    fractions (phi, cs, fs);
    fprintf(ferr, "STL file is saved in fs... \n");
    fclose (fp);
}

u.n[left] = dirichlet(v_in);
u.t[left] = dirichlet(0);
u.r[left] = dirichlet(0);
//p[left] = neumann(0);
//pf[left] = neumann(0);

u.n[right] = neumann(0);
u.t[right] = neumann(0);
p[right] = dirichlet(0);
pf[right] = dirichlet(0);
/**
The boundary condition is zero velocity on the embedded boundary. */

u.n[embed] = dirichlet(0);
u.t[embed] = dirichlet(0);
u.r[embed] = dirichlet(0);

int main(int argc, char * argv[])
{
    if (argc > 1) maxlevel = atoi(argv[1]); //convert from string to int
//    origin (-0.5, -0.5, -0.5);
    size(1.0);
    foreach_dimension()
    periodic (top);
    periodic (front);

    /**
    We turn off the advection term. The choice of the maximum timestep
    and of the tolerance on the Poisson and viscous solves is not
    trivial. This was adjusted by trial and error to minimize (possibly)
    splitting errors and optimize convergence speed. */

    stokes = true;
    DT = 1e-9;
    TOLERANCE = 1e-4;
    NITERMIN = 2;
    N = 1 << 7;
    v_in =
    const face vector muc[] = {MU,MU,MU};
    mu = muc;
    (const) scalar rhoc[] = RHO;
    rho = rhoc;
    run();
}

scalar un[];

event init (t = 0) {
    if (!restore (file = "restart")) {
        if (pid() > 0) {
            fprintf(ferr, "you have run in parallel reading ONLY in serial N parallel=%d.. "
                          "Run in serial save dump, then rename it to restart, then run in parallel", pid());
            exit(1992);
        }
        /**
        We define the bone embedded geometry. */

        bone(cs, fs);
        /**
        We initialize the reference velocity. */

        foreach() un[] = u.x[];
        DT = 1e-9;
    }
}

event set_dtmax (i++) {
    DT *= 1.05;
    DT = min(DT, CFL/pow(2, maxlevel+3));
    fprintf(ferr, "set_dtmax: tnext= %g DT= %g", tnext, DT);
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

/**
We check for a stationary solution. */

event logfile (i++; i <= 500)
{
    double avg = normf(u.x).avg, du = change (u.x, un)/(avg + SEPS);
    stats s = statsf (u.x);
    fprintf (ferr, "maxlevel= %d t= %g dt= %g i= %d avg= %g mgp.i= %d mgp.nrelax= %d mgp.minlevel= %d mgu.i= %d mgu.nrelax= %d mgu.minlevel= %d du= %.3g mgp.resa*dt= %.3g mgu.resa= %.3g u.x.sum()= %.3g normf(p).max= %.3g\n",
    maxlevel, t, dt, i, avg,
    mgp.i, mgp.nrelax, mgp.minlevel,
    mgu.i, mgu.nrelax, mgu.minlevel,
    du, mgp.resa*dt, mgu.resa, s.sum, normf(p).max);

    /**
    If the relative change of the velocity is small enough we stop this
    simulation. */

    if (i > 1 && (avg/dt < 1e-4 || du/dt < 0.01)) {
        /**
        We are interested in the permeability $k$ of the medium, which is
        defined by
        $$
        U = \frac{k}{\mu}\nabla p = \frac{k}{\mu}\rho g
        $$
        with $U$ the average fluid velocity. */
        fprintf (ferr, "converges... t= %g i= %d vtk_file= %d maxlevel= %d kappa= %g\n", t, i, iter_fp, maxlevel, s.sum/s.volume);
//        event ("vtk_file");
        boundary (all); // this is necessary since BCs depend on embedded fractions
        event ("adapt");
        event ("snapshot");
        fprintf(ferr, "Finish\n");
        return 0;
    }else{
        fprintf (ferr, "no convergence... t= %g i= %d vtk_file= %d maxlevel= %d kappa= %g\n", t, i, iter_fp, maxlevel, s.sum/s.volume);
    }
}

//event vtk_file (i += 50){
event vtk_file (t += 0.01){
    char subname[80]; sprintf(subname, "bone");
    scalar l[];
    foreach() {l[] = level;}
    output_vtu_MPI( (scalar *) {cs, p, l}, (vector *) {u}, subname, 1 );
}

event snapshot (i +=100) {
    char name[80];
    sprintf (name, "dump-%d", i);
    dump (file = name);
    fprintf(ferr,"dumped file=%s in t= %g dt= %g i= %d\n", name, t, dt, i);
}

/**
## Mesh adaptation

This computation is only feasible thanks to mesh adaptation, based
both on volume fraction and velocity accuracy. */
#define ADAPT_SCALARS {cs, u}
#define ADAPT_EPS_SCALARS {cseps, ueps, ueps, ueps}
event adapt (i++){
//    double eps_arr[] = ADAPT_EPS_SCALARS;
//    MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) {cs, u}, {cseps, ueps, ueps, ueps}, maxlevel = maxlevel, minlevel = minlevel);
}
/**
![Boundary of the bone medium. Fluid is flowing inside this
 volume.](bone3D/cs-8.png)

![Isosurface of the x-component of the velocity.](bone3D/nu-8.png)

![Cross-section at $z = 0$ showing the mesh and norm of the
 velocity.](bone3D/cross-8.png)

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
plot '../bone3D.ref' u 2:9 w l t '', '' u 2:10 w l t '', \
    '' u 2:11 w l t '', '' u 2:12 w l t '', '' u 2:13 w l t '', \
    'log' u 2:9 w p t 'du', '' u 2:10 w p t 'resp', \
    '' u 2:11 w p t 'resu', '' u 2:12 w p t 'u.x.sum', '' u 2:13 w p t 'p.max'
~~~

## See also

* [Stokes flow past a periodic array of spheres](/src/test/spheres.c)
* [Stokes flow past a periodic array of cylinders](/src/test/cylinders.c)
* [Stokes flow through a complex bone medium](/src/test/bone.c)
*/

/**
We output fields and dump the simulation. */

//        scalar nu[];
//        foreach()
//        nu[] = norm(u);
//        boundary ({nu});
//
//        view (fov = 32.2073, quat = {-0.309062,0.243301,0.0992085,0.914026},
//                tx = 0.0122768, ty = 0.0604286, bg = {1,1,1},
//                width = 600, height = 600);
//        box();
//        draw_vof("cs", "fs", fc = {0.5,0.5,0.5});
//        char name[80];
//        sprintf (name, "cs-%d.png", maxlevel);
//        save (name);
//
//        box();
//        isosurface ("u.x", 1e-5, fc = {0.,0.7,0.7});
//        sprintf (name, "nu-%d.png", maxlevel);
//        save (name);
//
//        view (fov = 19.1765, quat = {0,0,0,1},
//                tx = 0.0017678, ty = 0.00157844, bg = {1,1,1},
//                width = 600, height = 600);
//        squares ("nu", linear = true);
//        cells();
//        sprintf (name, "cross-%d.png", maxlevel);
//        save (name);
//
//        sprintf (name, "dump-%d", maxlevel);
//        dump (name);
