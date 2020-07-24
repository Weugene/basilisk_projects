#define BRINKMAN_PENALIZATION 1
#define DEBUG_MINMAXVALUES
//#define DEBUG_BRINKMAN_PENALIZATION 1
//#define DEBUG_OUTPUT_VTU_MPI
#define FILTERED
#define STOKES
scalar fs[];
scalar my_kappa[];
#include "grid/octree.h"
#include "../src_local/centered-weugene.h"
#include "two-phase.h"
#ifdef STOKES
    #include "navier-stokes/conserving.h"
#endif
#include "tension.h"
#include "../src_local/adapt_wavelet_limited.h"
#include "../src_local/utils-weugene.h"
#include "../src_local/output_vtu_foreach.h"
#include "view.h"
#include "maxruntime.h"


#define AIR_WATER
//#define AIR_GLYCEROL
#define uexact(x,y,z) 2.*(1. - 4*sq(y) - 4*sq(z))
//#define uexact(x,y,z) 0.25*(G/mu1)*(sq(0.5) - sq(y) - sq(z))
//Channel cross section Lyy*Lzz
double Vd, deq, dst = 0.2, rst = 0.1;
double RhoR, MuR;
#if defined(AIR_WATER)
    double Rho1 = 997, Rho2 = 1.204;
    double Mu1 = 0.88e-3, Mu2 = 0.019e-3;
    double Sigma = 72.8e-3;
    double diam_tube = 514e-6;
#elif defined(AIR_GLYCEROL)
    double Rho1 = 1250, Rho2 = 1.204;
    double Mu1 = 550e-3, Mu2 = 0.019e-3;
    double Sigma = 63.4e-3;
    double diam_tube = 494e-6;
#endif
double Ca; // Ca = Mu*Ud/sigma
double Ca_mod; // Ca_mod = Mu*Umean/sigma
double Re; //Reynolds
double G;
double Umean, vol_bubble, r_bub, l_bub;
double x_init = 2;
int maxlevel = 13;
int minlevel = 5;
int LEVEL = 6;
int adapt_method = 1;
double fseps = 1e-3, ueps = 1e-3;
scalar un[];

int main (int argc, char * argv[]) {
    maxruntime (&argc, argv);
    size(40);
    origin(0., -L0/2., -L0/2.);
    init_grid(1 << LEVEL);
    eta_s = 1e-5;
    TOLERANCE = 1e-6;
    NITERMIN = 1;
    NITERMAX = 100;
    DT = 1e-3;
    relative_residual_poisson = true;
    relative_residual_viscous = true;
// Case 9e
//    Ca = 0.163; //Ca = Ud*Mu1/sigma
//    Umean = 0.01145; // m/s
//    Vd = 0.0780e-9; // m^3

// Case 10e
    Ca = 0.023; //Ca = Ud*Mu1/sigma
    Umean = 1.580; // m/s
    Vd = 0.2179e-9; // m^3

    if (argc > 1)
        maxlevel = atoi (argv[1]);
    if (argc > 2)
        Ca = atof (argv[2]);
    if (argc > 3)
        Umean = atof (argv[3]);
    if (argc > 4)
        Vd = atof (argv[4]);


    deq = pow(6*Vd/pi, 1./3.);// 0.0005301091821 m
    dst = deq/diam_tube;// 1.0730955104
//    Umean = G*sq(0.5)/(8*Mu1);
    Ca_mod = Mu1*Umean/Sigma;
    Re = Umean*diam_tube*Rho1/Mu1;
    G = 32.0*Mu1*Umean/sq(diam_tube);
    rst = 0.5*dst;
    vol_bubble = (4./3.)*pi*cube(rst);
    r_bub = min(rst, 0.4);
    l_bub = cube(rst)/sq(r_bub);
    x_init = 1.7*l_bub;

    fprintf(ferr,"BP:             eta_s=%g,     TOLERANCE=%g    DT=%g\n"
                 "ADAPT:          minlevel=%d,  maxlevel=%d\n"
                 "Properties(SI): Mu1=%g Mu2=%g Rho1=%g Rho2=%g Sigma=%g G=%g  Umean=%g\n"
                 "Apparatus:      diam_tube=%g  length=%g\n"
                 "Bubble:         Vd=%g deq=%g\n",
                 eta_s, TOLERANCE, DT,
                 minlevel, maxlevel,
                 Mu1, Mu2, Rho1, Rho2, Sigma, G, Umean,
                 diam_tube, L0,
                 Vd, deq);
    // Dimensionless parameters:
    // Averaging on diam_tube=1 and Umean=1, Mu1=1 and Rho1=1 p' = p/(Rho1*Umean^2)
//    G /= Mu1*Umean/sq(diam_tube);
    G /= Rho1*sq(Umean)/diam_tube;
    Umean /= Umean;
    diam_tube /= diam_tube;
    RhoR = Rho1/Rho2;
    MuR = Mu1/Mu2;

    rho1 = 1.;// water
    rho2 = 1./RhoR; // air
    mu1 = 1./Re;
    mu2 = 1./(MuR*Re);
    f.sigma = 1./(Re*Ca_mod);
    fprintf(ferr,"Dimensionless Parameters: mu1=%g mu2=%g rho1=%g rho2=%g sigma=%g G=%g  Umean=%g\n"
                 "Dimensionless nums:       Re=%g  Ca=%g  Ca_mod=%g\n"
                 "Bubble:                   vol_bubble=%g dst=%g  rst=%g  r_bub=%g l_bub=%g\n",
            mu1, mu2, rho1, rho2, f.sigma, G, Umean,
            Re, Ca, Ca_mod,
            vol_bubble, dst, rst, r_bub, l_bub);
    run();
}

//BCs

//Inflow
u.n[left] = dirichlet(uexact(x,y,z)*(1-fs[]));
p[left] = neumann(0.);
pf[left] = neumann(0.);

//Outflow
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);
pf[right] = dirichlet(0.);

void geometry(scalar fs)
{
    //Channel
    vertex scalar phi[];
    foreach_vertex() {
        phi[] = sq(y) + sq(z) - sq(0.5);
    }
    boundary ({phi});
    fractions (phi, fs);
}

void bubble(scalar f)
{
    vertex scalar phi[];
    foreach_vertex() {
        phi[] = sq((x - x_init)/l_bub) + sq(y/r_bub) + sq(z/r_bub) - sq(1);
    }
    boundary ({phi});
    fractions (phi, f);
}

int maXlevel(double x,double y, double z){
    double x0 = fabs(x - x_init - Umean*t);
    int n = ceil(max(0, 2*(x0/l_bub - 2)));
    return max(maxlevel-n, 8);
}

event init (t = 0) {
  if (!restore (file = "restart")) {
    int it = 0;
    do {
      it++;
      geometry(fs);
      bubble(f);
    }while (adapt_wavelet_limited({fs, f}, (double []){fseps, fseps},
                                  maXlevel, minlevel).nf != 0 && it <= 10);

    foreach() {
      u.x[] = (1 - fs[])*uexact(x,y,z); //Init velocity
      u.y[] = 0;
      u.z[] = 0;
      un[] = u.x[];
      p[] = (1-f[])*2.0*f.sigma/rst;
      g.x[] = 0.5*((p[0] - p[-1])/Delta + (p[1] - p[0])/Delta);
    }
    event("vtk_file");
  }
}

//
//event logfile (i += 10; t <= 10) {
//  double du = change (u.x, un);
//  fprintf (ferr, "%d %g %g\n", i, t, du);
//  fflush (ferr);
//  if (i > 0 && du < 1e-6) //Convergence criteria
//    return 1; //Stop the simulation
//}

event logfile (i +=100)
{
    double avggas = L0*1 - normf(f).avg;
    scalar umag[];
    double xcg = 0, volume = 0, volumeg = 0, velamean = 0, velgx = 0, velgy = 0, dvtmp;
    foreach(reduction(+:xcg) reduction(+:volume) reduction(+:volumeg)
            reduction(+:velgx) reduction(+:velgy)
            reduction(+:velamean) ) {
        if (fs[]<1){
            dvtmp = (1.0 - f[])*(1.0 - fs[])*dv(); // gas volume
            volumeg += dvtmp;//gas liquid
            volume += (1.0 - fs[])*dv();//channel volume
            umag[] = norm(u); // the length of u
            xcg   += x*dvtmp;// Along x
            velamean += umag[]*dvtmp;//mean velocity of gas
            velgx += u.x[]*dvtmp;//mean velocity of gas Ox
            velgy += u.y[]*dvtmp;//mean velocity of gas Oy
        }
    }
    xcg /= volumeg; velgx /= volumeg; velgy /= volumeg; velamean /= volumeg;
    //norm statu = normf_weugene(umag, fs); // outputs avg, rms, max, volume
    fprintf (ferr, "maxlevel= %d i= %d t= %g dt= %g avggas= %g velgx= %g valgy= %g velamean= %g velgx/U0-1= %g xcg= %g\n",
            maxlevel, i, t, dt, avggas, velgx, velgy, velamean, (velgx/Umean - 1), xcg);
}

//event profiles (t = end) {
//  scalar * my_list = {u.x}; //list of scalars I want to export
//  int len_my_list = list_len(my_list); //lenght of my list
//  int np = 100;
//  double v[(np+1)*len_my_list]; //number of interpolated points (np+1) times number of scalars (len_my_list)
//
//  //line 1
//  coord a[np+1];
//  for (int n = 0; n <= np; n++) {
//    a[n].x = 20.;
//    a[n].y = 0.;
//    a[n].z = -0.5 + (1./np) * n;
//  }
//  interpolate_array(my_list, a, np+1, v, true);
//
//  if (pid()==0) {
//    FILE * fp1 = fopen ("profile_1", "w");
//    for (int n = 0; n <= np; n++) {
//      fprintf(fp1, "%g %g %g %g\n", a[n].x, a[n].y, a[n].z,
//	      v[n*len_my_list]);
//    }
//    fclose(fp1);
//  }
//}

//event movie (t += 0.1) {
//  view (fov = 22.4578, quat = {-0.707107,-0,-0,0.707107}, tx = -0.5, ty = 0.,
//  	bg = {0.3,0.4,0.6}, width = 600, height = 600, samples = 1);
//
//  clear();
//  squares("u.x", min=0, max=1.5, alpha = 0, n = {0,1,0});
//  save("movie_ux.mp4");
//}

event snapshot (i += 1000)
{
  char name[80];
  sprintf(name, "dump-%04g",t);
  dump (file = name);
}

void exact(vector ue)
{
    foreach() {
        ue.x[] = (1 - fs[])*uexact(x,y,z);
        ue.y[] = 0;
        ue.z[] = 0;
    }
    boundary((vector *){ue});
}
event vtk_file (t += 0.1)
{
  char subname[80]; sprintf(subname, "tube_bp");
  vector ue[]; exact(ue);
  scalar l[]; foreach() l[] = level;
//  vector du[]; foreach() foreach_dimension() du.x[] = u.x[] - ue.x[]; boundary((vector *){du});
  output_vtu_MPI( subname, t + dt, (scalar *) {p, fs, f, my_kappa, l}, (vector *) {u, ue} );
}



#define ADAPT_SCALARS {fs, f, u}
#define ADAPT_EPS_SCALARS {fseps, ueps, ueps, ueps}
event adapt (i++)
{
    double eps_arr[] = ADAPT_EPS_SCALARS;
    MinMaxValues (ADAPT_SCALARS, eps_arr);
    if (adapt_method)
        adapt_wavelet_limited  ((scalar *) {fs, f, u}, (double []){fseps, fseps, ueps, ueps, ueps}, maXlevel, minlevel);
    else
        adapt_wavelet ((scalar *) {fs, f, u}, (double []){fseps, fseps, ueps, ueps, ueps}, maxlevel = maxlevel, minlevel = minlevel);
    geometry(fs);
}

event stop(t=10);