/**
#  Encapsulation of water drop by Silicon oil pool

<p align="center">
<video width="20%" controls>
  <source src="https://www.dropbox.com/s/59tlyvx8rppp7f0/Experiment.mp4?dl=1" type="video/mp4">
  <caption><p align="center">Experiment: As the water drop reaches the liquid (oil)
   bath, it attains a quasi-static shape with the help of a thin air-layer.
   When the air drains out, a triple contact line (between air, water, and oil) forms which generates capillary waves and destabilizes the system.
   Owing to the positive spreading coefficient of the oil, it encapsulates the water droplet.
   Experiment done by [Utakarsh Jain](https://pof.tnw.utwente.nl/people/profile/1169)
   and [Mazi (Maziyar) Jalaal](https://pof.tnw.utwente.nl/people/profile/1169).</caption>
</video>

*/


// 1 is Si Pool, 2 is Water Drop and 3 is air
#include "../src_local/centered-weugene.h"
// Smering out the interfaces for high density and viscosity ratios.
//#define FILTERED
#include "../src_local/three-phase-weugene.h"
#include "../src_local/surface-tension-decompose.h"
#include "distance.h"

#define MAXlevel 9                                             // maximum level
#define MINlevel 4                                              // maximum level
#define tmax 15.0                                                // maximum time
#define tsnap (1e-2)
// Error tolerances
#define fErr (1e-3)                                  // error tolerance in VOF
#define K1Err (1e-1)                                 // error tolerance in KAPPA 1e-3
#define K2Err (1e-1)                                 // error tolerance in KAPPA
#define VelErr (2e-1)                                // error tolerances in velocity
#define OmegaErr (1e-2)                              // error tolerances in velocity
#define EPS_MAXA 1                                   // method of eps calculation

/**
Different scales used in the problem: $U_\sigma \equiv \sqrt{{\sigma}/{\rho_l R_0}}$
and time: $t_\sigma \equiv {R_0}/{U_\sigma}$.<br/>
Navier-Stokes equation is therefore:
$$ \hat{\rho}\left(\frac{\partial U_i}{\partial t} + U_j\frac{\partial U_i}{\partial X_j}\right) = -\frac{\partial P}{\partial X_i} + \frac{\partial}{\partial X_j}\left(2\left(\hat{\mu}Oh\right)D_{ij}\right) + \kappa\delta_s\hat{n}_{i} + Bo\hat{g}_i $$
*/
#define La (3.6e5)
#define Oh 1
//#define Oh sqrt(1./La)
#define Bo 0
#define Mu12 1
#define Mu32 1
#define Rho12 1.0
#define Rho32 1


/** Surface tesnions
Non-dimensioalized using the surface tenison coefficient of $\sigma_{23}$.<br/>
Notice that $\sigma_{12}/\sigma_{23} + \sigma_{13}/\sigma_{23} < 1$. The spreading coefficient
of liquid pool (1) is positive. It will try to maximize its surface area.
\cos \theta = \frac{\sigma_{23}-\sigma_{12}}{\sigma_{13}} =\frac{\sigma_{sg}-\sigma_{ls}}{\sigma_{lg}}
 (1-1.7)/1=-0.7
 */
#define SIGMA12by23 (0.1)
#define SIGMA13by23 (1)

// density
#define Rho1 Rho12                                         // density of phase 1
#define Rho2 (1.0)                                       // density of phase 2
#define Rho3 Rho32                                         // density of phase 3

// viscosity
#define Mu1 (Mu12*Oh)                            // Dynamic viscosity of phase 1
#define Mu2 (1.00*Oh)                            // Dynamic viscosity of phase 2
#define Mu3 (Mu32*Oh)                            // Dynamic viscosity of phase 3

// domain
#define Ldomain 1                                 // Dimension of the domain

// boundary conditions
u.n[right] = neumann(0.);
u.n[left] = neumann(0.);
u.n[top] = neumann(0.);
u.n[bottom] = dirichlet(0.);

f[left] =neumann(0);
f[right]=dirichlet(0);
f[top] =dirichlet(0);
f[bottom] =dirichlet(0);

fs[left] =neumann(0);
fs[right]=neumann(0);
fs[top] =dirichlet(0);
fs[bottom] =dirichlet(1);
//p[left] = dirichlet(0.);
//p[right] = dirichlet(0.);
p[top] = dirichlet(0.);
//p[bottom] = dirichlet(0.);
scalar omega[];
int main()
{
    init_grid (1024);
    rho1 = Rho1; mu1 = Mu1;
    rho2 = Rho2; mu2 = Mu2;
    rho3 = Rho3; mu3 = Mu3;
    sigma13 = SIGMA13by23;
    sigma23 = 1.0;
    sigma12 = SIGMA12by23;

    L0=Ldomain;
    X0=0;
    Y0=-L0/8;
    stokes = true; //added
    DT = 1e-3;
#if TREE
    for (scalar s in {f, fs}){
    s.refine = s.prolongation = fraction_refine;
    boundary ({s});
  }
#endif
    run();
}
event init(t = 0) {
    if (!restore (file = "restart")) {
        scalar tmp[];
        vector gtmp[];
        int iter = 0;
        do {
            iter++;
            foreach() {
                fs[] = (y < 0) ? 1 : 0;
                f[] = (sq(x) + sq(y+L0/4) - sq(L0/2) < 0) && (y > 0) ? 1 : 0;
                tmp[] = f[] + 2*fs[];
                u.x[] = 0.0;
                u.y[] = 0.0;
            }
            vorticity (u, omega);
            boundary ((scalar *){omega});
            boundary ({f, fs, u, tmp});
            gradients ({tmp}, {gtmp});
        }while (adapt_wavelet({gtmp.x, gtmp.x, u.x, u.y, omega}, (double []){fErr, fErr, VelErr, VelErr, OmegaErr},
                              maxlevel = MAXlevel, minlevel=MINlevel).nf != 0 && iter <= 15);

        fprintf(stderr, "init refinement iter=%d", iter);
    }else{
        fprintf(stderr, "RESTART from file");
    }

    event("properties");
    event("acceleration");
    event ("vtk_file");
}
// Gravity
event acceleration(i++){
face vector av = a;
foreach_face(y){
        av.y[] -= Bo;
}
}

event velocity_correction(i++){
    foreach() foreach_dimension() u.x[] *= (1.0 - fs[]);
    boundary({u});
}
//Output
#include "../src_local/output_vtu_foreach.h"
event vtk_file (i += 100)
{
    char subname[80]; sprintf(subname, "enc");
    scalar l[]; foreach() l[] = level;

    vorticity (u, omega);
    output_vtu_MPI( (scalar *) {l, f, fs, rho, p, pf, phi1, phi2, omega}, (vector *) {a, u, uf, mu}, subname);
}

#define ADAPT_SCALARS {ga.x, ga.y}
#define ADAPT_EPS_SCALARS {fErr, fErr}

event adapt(i++)
{
/**
We use AMR based on the curvature values of VOF fields, vorticity, and velocities.
*/
    vorticity (u, omega);
    scalar a[];
    vector ga[];
    foreach() a[] = f[] + 2*fs[];
    boundary({a});
    gradients ({a}, {ga});
    double eps_arr[] = ADAPT_EPS_SCALARS;
//    MinMaxValues(ADAPT_SCALARS, eps_arr);
    adapt_wavelet ((scalar *) ADAPT_SCALARS,
                    eps_arr, maxlevel = MAXlevel, minlevel = MINlevel);
}

event stop (t = tmax);
// Outputs
//event writingFiles (t = 0; t += tsnap; t <= tmax + tsnap) {
//  dump (file = "dump");
//  char nameOut[80];
//  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
//  dump (file = nameOut);
//}

//event logWriting (i++) {
//  double ke = 0.;
//  foreach (reduction(+:ke)){
//    ke += sq(Delta)*(sq(u.x[]) + sq(u.y[]))*rho(f[],fs[]);
//  }
//  static FILE * fp;
//  if (i == 0) {
//    fprintf (ferr, "i dt t ke\n");
//    fp = fopen ("log", "w");
//    fprintf (fp, "i dt t ke\n");
//    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
//    fclose(fp);
//  } else {
//    fp = fopen ("log", "a");
//    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
//    fclose(fp);
//  }
//  fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);
//}

/**
## Running the code

Use the following `run.sh` script

~~~bash
#!/bin/bash
qcc -O2 -Wall Encapsulation.c -o Encapsulation -lm
./Encapsulation
~~~

# Output and Results
The post-processing codes and simulation data are available at: [PostProcess](https://www.dropbox.com/sh/1lvcuf20vcuj8xw/AACRDVOu8nM-a-gYpqbMgarQa?dl=0)

## The process
<p align="center">
<video width="50%" controls>
  <source src="https://www.dropbox.com/s/m0g5ql45uiviqws/BasiliskWebsite.mp4?dl=1" type="video/mp4">
  <caption><p align="center">Temporal evolution of the droplet encapsulation process.
  The triple-contact line generates a train of capillary waves which travel through the droplet interface
  and helps in the encapsulation process. In each image, the left part shows vorticity, and the right shows velocity magnitude.</caption>
</video>
## Facets to show the triple point
<p align="center">
<video width="50%" controls>
  <source src="https://www.dropbox.com/s/0v6920zj4yb1706/TriplePoint.mp4?dl=1" type="video/mp4">
</video>
## AMR
<p align="center">
<video width="50%" controls>
  <source src="https://www.dropbox.com/s/co9sm7vdwdgb5a6/Bview.mp4?dl=1" type="video/mp4">
  <caption><p align="center">The maximum refinement for the simulation is 10. This means that $R_0/\Delta$ is 128. Left size of this video show the evolution of the adaptive mesh refinement grids. The colormap on the left shows different mesh resolution. On the right, vorticity is shown. (Fix me! Need a better method to appreciate AMR!)</caption>
</video>

# Bibliography

* [Duchemin, Laurent. Quelques problemes fortement non-linéaires de surface libre et leur résolution numérique. PhD theis, Université de la Méditerranée-Aix-Marseille II, (2001).](https://tel.archives-ouvertes.fr/tel-00084132/document)

* Wong, Clint YH, Mokhtar Adda-Bedia, and Dominic Vella. Non-wetting drops at liquid interfaces: From liquid marbles to Leidenfrost drops. Soft matter 13, no. 31 (2017): 5250-5260.
[doi: 10.1039/C7SM00990A](https://pubs.rsc.org/en/content/articlelanding/2017/sm/c7sm00990a#!divAbstract)
*/
