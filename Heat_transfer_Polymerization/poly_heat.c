#define BRINKMAN_PENALIZATION 1
#define DEBUG_MINMAXVALUES
//#define DEBUG_BRINKMAN_PENALIZATION
#define DEBUG_MODE_POISSON
//#define DEBUG_OUTPUT_VTU_MPI
#define FILTERED
#define JACOBI 1

//#define PRINT_ALL_VALUES
//#define STOKES

scalar omega[];
scalar l2[];
#include "../src_local/centered-weugene.h"
#include "../src_local/rheology_model.h"
#include "tension.h"
#include "../src_local/output_vtu_foreach.h"

double Uin, Tin;
double RhoR, RhoRS, MuR, MuRS, CpR, CpRS, KappaR, KappaRS;
double Rho1, Rho2, Rho3;
double Mu1, Mu2, Mu3;
double Kappa1, Kappa2, Kappa3;
double CP1, CP2, CP3;
double Sigma;
double cyl_diam, domain_size;
int Ncx, Ncy; //number of cylinders along Ox, Oy
double Ca; // Ca = Mu*Ud/sigma
double Re; //Reynolds
double G;
double Umean;
double x_init = 2;
int maxlevel = 10;
int minlevel = 5;
int LEVEL = 9;
int adapt_method = 1; // 0 - traditional, 1 - using limitation, 2 - using array for maxlevel
int snapshot_i;
double dt_vtk;
double feps=1e-3, fseps = 1e-3, ueps = 1e-2, Teps = 1e-3, aeps = 1e-3;
double TOLERANCE_P = 1e-5, TOLERANCE_V = 1e-5;

int main(int argc, char * argv[]) {
    L0 = 8.;
    origin (-L0/2, -L0/2.);
	eta_s = 1e-5;
    TOLERANCE = 1e-6;
    NITERMIN = 1;
    NITERMAX = 100;
    N = 1 << LEVEL;
    CFL = 0.1;
    DT = 1e-5;
	dt_vtk = 1e-20;
	snapshot_i = 100;
    periodic(top);
	if (argc > 1)
        maxlevel = atoi (argv[1]);
	if (argc > 2)
        iter_fp = atoi (argv[2]);
    if (argc > 3)
        dt_vtk = atof (argv[3]);
    if (argc > 4)
        snapshot_i = atoi (argv[4]);
    stokes = true;
// Physical parameters
	Uin = 1e-2, Tin = 300;
	Ncx = 5, Ncy = 5;
	cyl_diam = 20e-6;//m 5-25 microns
    Rho1 = 1140, Rho2 = 1, Rho3 = 2000;
    Mu1 = 0.155, Mu2 = 1.81e-5, Mu3 = 1;
	CP1 = 1100, CP2 = 1006, CP3 = 840;//J/(kg*K)
	Kappa1 = 1100, Kappa2 = 0.02535, Kappa3 = 1.11;//W/(m*K)
	Htr = 1;
	Arrhenius_const = 4.53e+7;//1/s
	Ea_by_R = 72900/8.314;// Kelvin
	Arrhenius_const = 10;//1/s
	n_degree = 1.667;
	m_degree = 0.333;
	Sigma = 0.072;// N/m
 fprintf(ferr,   "BP:        eta_s=%g,     DT=%g\n"
                 "Solver:    NITERMIN=%d,  NITERMAX=%d, TOLERANCE=%g  relative_residual_poisson=%d relative_residual_viscous=%d\n"
                 "OUTPUT:    dt_vtk=%g.    number of procs=%d\n"
                 "ADAPT:     minlevel=%d,  maxlevel=%d, feps=%g, fseps=%g, ueps=%g, Teps=%g, aeps=%g\n"
                 "Props(SI): Mu1=%g, Mu2=%g, Mu3=%g, Rho1=%g, Rho2=%g,  Rho3=%g\n,"
                 "           Kappa1=%g, Kappa2=%g, Kappa3=%g, CP1=%g, CP2=%g, CP3=%g\n"
				 "			 Sigma=%g Uin=%g Tin=%g Htr=%g Arrenius=%g Ea_by_R=%g\n"
                 "Apparatus: cyl_diam=%g  domainSize=%g\n",
                 eta_s, DT,
                 NITERMIN, NITERMAX, TOLERANCE, relative_residual_poisson, relative_residual_viscous,
                 dt_vtk, npe(),
                 minlevel, maxlevel, feps, fseps, ueps, Teps, aeps,
                 Mu1, Mu2, Mu3, Rho1, Rho2, Rho3, 
                 Kappa1, Kappa2, Kappa3, CP1, CP2, CP3,
                 Sigma, Uin, Tin, Htr, Arrhenius_const, Ea_by_R,
                 cyl_diam, domain_size);
// Dimensionless numbers
	Re = Uin*cyl_diam*Rho1/Mu1;
	Ca = Mu1*Uin/Sigma;
// Dimensionless parameters are chosen cyl_diam, rho1, Cp1, Tin, Uin
	RhoR = Rho2/Rho1, RhoRS = Rho3/Rho1;
	MuR = Mu2/Mu1, MuRS = Mu3/Mu1;
	CpR = CP2/CP1, CpRS = CP3/CP1;
	KappaR = Kappa2/Kappa1, KappaRS = Kappa3/Kappa1;
    rho1 = 1; rho2 = RhoR; rho3 = RhoRS;
    mu1 = 1/Re; mu2 = mu1*MuR; mu3 = mu1*MuRS;
	Cp1 = 1, Cp2 = CpR, Cp3 = CpRS;//J/(kg*K)
	kappa1 = Kappa1/(Rho1*CP1*cyl_diam*Uin), kappa2 = kappa1*KappaR, kappa3 = kappa1*KappaRS;//W/(m*K)
	Htr /= CP1*Tin;
	Arrhenius_const *= cyl_diam/Uin;
	Ea_by_R /= Tin;
    f.sigma = 1./(Re*Ca);
	
	Uin = 1;
	Tin = 1;
	fprintf(ferr,"Dim-less vars: mu1=%g mu2=%g mu3=%g rho1=%g rho2=%g rho3=%g\n" 
				 "			 	 kappa1=%g kappa2=%g kappa3=%g Cp1=%g Cp2=%g Cp3=%g\n"
				 "               sigma=%g  Uin=%g Tin=%g\n"
                 "Dim-less nums: Re=%g  Ca=%g\n",
				 mu1, mu2, mu3, rho1, rho2, rho3, 
				 kappa1, kappa2, kappa3, Cp1, Cp2, Cp3,
				 f.sigma, Uin, Tin,
                 Re, Ca);
    run();
}
//#define T_BC (0.5*(TMAX + TMIN) + 0.5*(TMAX - TMIN)*tanh((x)/(L0/10)))

u.n[left]  = dirichlet(Uin);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(1);
T[left]    = dirichlet(Tin);
fs[left]   = dirichlet(0);
alpha_doc[left] = neumann(0);//inflow is fresh resin

u.n[right] = neumann(0);
p[right]   = dirichlet(0);
pf[right]  = dirichlet(0);
T[right]    = neumann(0);
fs[right]   = neumann(0);
alpha_doc[right] = neumann(0);


/*double cubeF(double x, double y, double z, coord center, double size) {
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
            foreach_dimension()
            pnt_dist.x = centers[i].x - centers[j].x;
            if (mynorm(pnt_dist) < dist) {
                i--;
                break;
            };
        }
    }



    double phi = HUGE;
    for (int i = 0; i < ns; i++) {
        printf("i=%d x=%g y=%g R=%g\n", i, centers[i].x, centers[i].y, R[i] );
        foreach_dimension()
            pnt_dist.x = mypoint.x - centers[i].x;
        phi = min(phi, (mynorm(pnt_dist) - R[i]));
    }
    return min(phi, -L0/8.0 - x);// with front
//    return phi; // no front
}

coord center={0,0,0};
//const int Ncyl=7;
//double size_box;
//double R;//0.5
//double dist;

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

double geometry(double x, double y, double z, const int Nlayer) {
    int icyl;
    double dist = 0.5*L0/Nlayer;
    double geom = cylindersOy(x, y, z, 0.0, Ncyl);

    for (icyl=1; icyl<Nlayer; icyl++){
        geom = max(cylindersOy(x, y+dist/2.*(icyl%2), z, icyl*dist, Ncyl + icyl%2), geom);
    }
    double c = cubeF(x, y, z, center, size_box);
    return min(geom,c);
}


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
        if (ellipse_shape)
            phi[] = sq((x - x_init)/l_bub) + sq(y/r_bub) + sq(z/r_bub) - sq(1);
        else if (cylinder_shape)
        {
            phi[] = intersection( 0.5*l_bub - fabs(x - x_init),  sq(r_bub) - sq(y) - sq(z) );
            phi[] = union(phi[], sq(r_bub) - sq(fabs(x - x_init) - 0.5*l_bub) - sq(y) - sq(z));
            phi[] *= -1;
        }
    }
    boundary ({phi});
    fractions (phi, f);
}
*/

event init (t = 0) {
    if (!restore (file = "restart")) {
        int iter = 0;
        do {
            iter++;
            foreach()
            {
                f[] = (sq(x+3) + sq(y) - sq(0.25) > 0 &&
                       sq(x+3.2) + sq(y-2) - sq(0.25) > 0 &&
                       sq(x+3) + sq(y-1) - sq(0.3) > 0 &&
                       sq(x+2.6) + sq(y+1.5) - sq(0.3) > 0 &&
                       sq(x+2.8) + sq(y+3) - sq(0.4) > 0) && (x < 0) ? 1 : 0;

            	fs[] = 0.5*(1+tanh((1 - sq(x) - sq(y))/0.01));
            }
            boundary ({f, fs});
        }while ((iter <=5) || ((adapt_wavelet({f, fs}, (double []){feps, fseps},
                              maxlevel = maxlevel, minlevel=minlevel).nf != 0) && (iter <= 15)));
        fprintf(stderr, "init refinement iter=%d", iter);
        foreach() {
            T[] = 1;
            alpha_doc[] = 0;
            u.x[] = Uin*(1-fs[]);
        }
        foreach_face(){
            kappav.x[] = var_hom(f[], fs[], kappa1, kappa2, kappa3);
        }
		boundary ({T, alpha_doc, u, kappav});
    }
	event("vtk_file");
}

//event adapt_step(i<=5)  DT = 1e-9;//event adapt_step(i<=5)  DT = 1e-9;




//event vtk_file (i += 1)
event vtk_file (t += dt_vtk)
{
	char subname[80]; sprintf(subname, "heat_pol");
    scalar l[]; foreach() l[] = level;
	output_vtu_MPI( subname, (iter_fp) ? t + dt : 0, (scalar *) {T, alpha_doc, p, fs, f, l}, (vector *) {u, a});
}


#if DUMP
event snapshot (i += snapshot_i)
//event snapshot (t += 1e-1)
{
    char name[80];
    sprintf(name, "dump-%04g",t);
    vorticity (u, omega);
    p.nodump = false;
    dump (file = name);
}
#endif

/**
We adapt according to the error on the embedded geometry, velocity and
tracer fields. */

event adapt (i++) {
	adapt_wavelet ({f, T, alpha_doc, u}, (double[]){feps, Teps, aeps, ueps, ueps}, maxlevel, minlevel);
}

event stop(t=L0/Uin);
/**
## See also

* [Same example with Gerris](http://gerris.dalembert.upmc.fr/gerris/examples/examples/cylinder.html)
*/
