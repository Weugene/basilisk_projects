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

int maxlevel = 8;
int minlevel = 5;
double cseps = 1e-3, ueps = 1e-3;
double RE, FROUDE;
//#define flux  1e-4 // 4-6 л\мин или 67-100 см^3\сек
#define S 4e-6 //cross section
#define RHO  1050 //1,050—1,060 g/cm³
#define MU 5e-3 //4-5 mPa×s
#define GRAVITY 9.8
#define v_in 0.1 // m/s in aort=0.5 m/s
#define L0_phys (sqrt(S))

//const face vector muv[] = {MU,MU,MU};
face vector muv[];
face vector alphav[];
scalar rhov[];

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

u.n[left] = dirichlet(v_in*cs[]);
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
    size(1);
    foreach_dimension()
        periodic (right);
    init_grid (32);
//    origin (-0.5, -0.5, -0.5);
//    origin (-0.2 -0.2 -0.2);


    /**
    We turn off the advection term. The choice of the maximum timestep
    and of the tolerance on the Poisson and viscous solves is not
    trivial. This was adjusted by trial and error to minimize (possibly)
    splitting errors and optimize convergence speed. */

    stokes = true;
    DT = 1e-3;
    TOLERANCE = 1e-4;
    NITERMIN = 1;
    RE = v_in*L0_phys*RHO/MU;
    FROUDE = v_in/sqrt(GRAVITY*L0_phys);
    fprintf(ferr, "Re=%g Fr=%g\n", RE, FROUDE);
    run();
}

scalar un[];

event defaults (i = 0) {
    const face vector grav[] = {1./sq(FROUDE), 0, 0};
    a = grav;
    mu = muv;
    alpha = alphav;
    rho = rhov;
}
event init (t = 0)
{
    if (!restore (file = "restart")) {
#ifdef _MPI
        fprintf(ferr, "you have run in parallel reading ONLY in serial N parallel=%d.. "
                      "Run in serial save dump, then rename it to restart, then run in parallel", pid());
        exit(1992);
#endif
        /**
        We define the bone embedded geometry. */
//        porous(cs, fs);
        bone(cs, fs);
        /**
        We initialize the reference velocity. */
        foreach() un[] = u.x[];
        DT = 1e-4;
    }
}
event properties (i++)
{
    foreach_face() {
        alphav.x[] = fm.x[]/RHO;
            face vector muv = mu;
            muv.x[] = fm.x[]*MU;
    }
    foreach()
        rhov[] = cm[]*RHO;

}

/**
We check for a stationary solution. */

event logfile (i++; i <= 500)
{
    double avg = normf (u.x).avg, du = change (u.x, un)/(avg + SEPS);
    stats s = statsf (u.x);
    fprintf (ferr, "maxlevel= %d t= %g dt= %g i= %d avg= %g mgp.i= %d mgp.nrelax= %d mgp.minlevel= %d mgu.i= %d mgu.nrelax= %d mgu.minlevel= %d du= %.3g mgp.resa*dt= %.3g mgu.resa= %.3g u.x.sum()= %.3g normf(p).max= %.3g\n",
    maxlevel, t, dt, i, avg,
    mgp.i, mgp.nrelax, mgp.minlevel,
    mgu.i, mgu.nrelax, mgu.minlevel,
    du, mgp.resa*dt, mgu.resa, s.sum, normf(p).max);

    /**
    If the relative change of the velocity is small enough we stop this
    simulation. */

    if (i > 10 && (avg/dt < 1e-4 || du/dt < 1e-2)) {
        /**
        We are interested in the permeability $k$ of the medium, which is
        defined by
        $$
        U = \frac{k}{\mu}\nabla p = \frac{k}{\mu}\rho g
        $$
        with $U$ the average fluid velocity. */
        fprintf (ferr, "converges... t= %g i= %d vtk_file= %d maxlevel= %d Uaver= %g kappa= %g avg/dt=%g du/dt=%g\n", t, i, iter_fp, maxlevel, s.sum/s.volume, (s.sum/s.volume)*MU/(RHO*GRAVITY), avg/dt, du/dt);
        event ("vtk_file");
        boundary (all); // this is necessary since BCs depend on embedded fractions
//        event ("adapt");
        event ("snapshot");
        fprintf(ferr, "Finish\n");
        return 0;
    }else{
        fprintf (ferr, "no convergence... t= %g i= %d vtk_file= %d maxlevel= %d kappa= %g\n", t, i, iter_fp, maxlevel, s.sum/s.volume);
    }
}

//event end_timestep (i += 50){
event end_timestep (t += 0.01)
{
    char subname[80]; sprintf (subname, "bone");
    scalar l[];
    foreach() {l[] = level;}
    output_vtu_MPI (subname, t+dt, (scalar *) {cs, p, l}, (vector *) {u} );
}

event snapshot (i +=100)
{
    char name[80];
    sprintf (name, "dump-%d", i);
    dump (file = name);
    fprintf (ferr,"dumped file=%s in t= %g dt= %g i= %d\n", name, t, dt, i);
}

/**
## Mesh adaptation

This computation is only feasible thanks to mesh adaptation, based
both on volume fraction and velocity accuracy. */
//#define ADAPT_SCALARS {cs, u}
//#define ADAPT_EPS_SCALARS {cseps, ueps, ueps, ueps}
//event adapt (i++)
//{
//    double eps_arr[] = ADAPT_EPS_SCALARS;
//    MinMaxValues (ADAPT_SCALARS, eps_arr);
//    adapt_wavelet ((scalar *) {u}, (double []){cseps, ueps, ueps, ueps}, maxlevel = maxlevel, minlevel = minlevel);
//}

event stop (t = 17);