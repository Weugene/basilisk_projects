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
//#include "axi.h"
#include "../src_local/centered-weugene.h"
// Smering out the interfaces for high density and viscosity ratios.
//#define FILTERED
#include "../src_local/three-phase-weugene.h"
#include "../src_local/tension_three-phase-weugene.h"
#include "distance.h"

#define MAXlevel 10                                             // maximum level
#define MINlevel 4                                              // maximum level
#define tmax 5.0                                                // maximum time
#define tsnap (1e-2)
// Error tolerances
#define fErr (1e-2)                                  // error tolerance in VOF
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
#define La (7.2e4)
#define Oh sqrt(1./La)
#define Bo 1.32
#define Mu12 20.0
#define Rho12 0.9630
#define Rho32 0.0010
#define Mu32  0.0185

/** Surface tesnions
Non-dimensioalized using the surface tenison coefficient of $\sigma_{23}$.<br/>
Notice that $\sigma_{12}/\sigma_{23} + \sigma_{13}/\sigma_{23} < 1$. The spreading coefficient
of liquid pool (1) is positive. It will try to maximize its surface area.
*/
#define SIGMA12by23 (1.0)
#define SIGMA13by23 (0.30)

// density
#define Rho1 Rho12                                         // density of phase 1
#define Rho2 (0.963)                                       // density of phase 2
#define Rho3 Rho32                                         // density of phase 3

// viscosity
#define Mu1 (Mu12*Oh)                            // Dynamic viscosity of phase 1
#define Mu2 (1.00*Oh)                            // Dynamic viscosity of phase 2
#define Mu3 (Mu32*Oh)                            // Dynamic viscosity of phase 3

// domain
#define Ldomain 8                                 // Dimension of the domain

// boundary conditions
u.n[right] = neumann(0.);
u.n[left] = neumann(0.);
u.n[top] = neumann(0.);
u.n[bottom] = dirichlet(0.);

//p[left] = dirichlet(0.);
//p[right] = dirichlet(0.);
//p[top] = dirichlet(0.);
p[bottom] = dirichlet(0.);
scalar KAPPA1[], KAPPA2[], omega[];
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
  X0=-L0/2;
  Y0=-L0/8;
  stokes = true; //added
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
        int iter = 0;
        do {
            iter++;
            foreach() {
                fs[] = (y < 0) ? 1 : 0;
                f[] = (sq(x) + sq(y) - sq(L0/8) < 0) && (y > 0) ? 1 : 0;
                f[] += fs[];
                u.x[] = 0.0;
                u.y[] = 0.0;
            }
            vorticity (u, omega);
            boundary ((scalar *){omega});
            boundary ({f, fs, u});
        }while (adapt_wavelet({f, fs, u.x, u.y, omega}, (double []){fErr, fErr, VelErr, VelErr, OmegaErr},
                              maxlevel = MAXlevel, minlevel=MINlevel).nf != 0 && iter <= 15);

        fprintf(stderr, "init refinement iter=%d", iter);
        foreach() {

        }
    }else{
        fprintf(stderr, "RESTART from file");
    }

    event ("vtk_file");
}
// Gravity
event acceleration(i++){
  face vector av = a;
  foreach_face(x){
    av.x[] -= Bo;
  }
}

event velocity_correction(i++){
    foreach() foreach_dimension() u.x[] *= (1.0 - fs[]);
    boundary({u});
}
//Output
#include "../src_local/output_vtu_foreach.h"
event vtk_file (i += 1)
{
    char subname[80]; sprintf(subname, "enc");
    scalar l[]; foreach() l[] = level;

    vorticity (u, omega);
    curvature(f, KAPPA1);
    curvature(fs, KAPPA2);

    output_vtu_MPI( (scalar *) {l, f, fs, rho, p, pf, KAPPA1, KAPPA2, omega}, (vector *) {u, uf, mu}, subname);
}

void MinMaxValues(const scalar * list, double * arr_eps) {// for each scalar min and max
    double arr[10][2];
    int ilist = 0;
    for (scalar s in list) {
        int mina= HUGE, maxa= -HUGE;
        foreach( reduction(min:mina) reduction(max:maxa) ){
            if (fabs(s[]) < mina) mina = fabs(s[]);
            if (fabs(s[]) > maxa) maxa = fabs(s[]);
        }
        arr[ilist][0] = mina;
        arr[ilist][1] = maxa;
        ilist++;
//        fprintf(stderr, "arr for i=%d", ilist);
    }

    for (int i = 0; i < ilist; i++){
#if EPS_MAXA == 1
        arr_eps[i] *=arr[i][1];
#else
        arr_eps[i] *= 0.5*(arr[i][0] + arr[i][1]);
#endif
        fprintf(stderr, "MinMaxValues: i=%d, min=%g, max=%g, eps=%g\n", i, arr[i][0], arr[i][1], arr_eps[i]);
    }


}

#define ADAPT_SCALARS {f, fs, omega}
#define ADAPT_EPS_SCALARS {fErr, 0.0001, OmegaErr}

event adapt(i++)
{
  /**
  We use AMR based on the curvature values of VOF fields, vorticity, and velocities.
  */
  vorticity (u, omega);
//  curvature(f, KAPPA1);
//  curvature(fs, KAPPA2);
//  foreach(){
////      KAPPA1[] = exp(-fabs(KAPPA1[]));
////      KAPPA2[] = exp(-fabs(KAPPA2[]));
////      KAPPA1[] = clamp(KAPPA1[], -1000, 1000);
////      KAPPA2[] = clamp(KAPPA2[], -1000, 1000);
//      KAPPA1[] = (fabs(KAPPA1[])>1000)? 1000: KAPPA1[];
////      KAPPA2[] = (fabs(KAPPA2[])>1000)? 0: KAPPA2[];
//  }
  double eps_arr[] = ADAPT_EPS_SCALARS;
  MinMaxValues(ADAPT_SCALARS, eps_arr);
  adapt_wavelet ((scalar *) ADAPT_SCALARS,
      eps_arr,
      maxlevel = MAXlevel, minlevel = MINlevel);

  fs.refine = fs.prolongation = fraction_refine;
//  boundary ({fs});


}
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
event stop (t = tmax);
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
