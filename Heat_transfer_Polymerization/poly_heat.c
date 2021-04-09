#define BRINKMAN_PENALIZATION 1
//#define DEBUG_MINMAXVALUES
//#define DEBUG_BRINKMAN_PENALIZATION
#define DEBUG_MODE_POISSON
#define REACTION_MODEL REACTION_MODEL_NON_AUTOCATALYTIC
//#define DEBUG_OUTPUT_VTU_MPI
#define FILTERED
//#define JACOBI 1
//#define RELATIVE_RES
//#define PRINT_ALL_VALUES
//#define STOKES

scalar omega[];
scalar l2[];
face vector av[];
scalar smoothed_f[], smoothed_fs[];
vector h[], hs[];
/**
# Contact angles

This file is used to impose contact angles on boundaries for
interfaces described using a [VOF](vof.h) tracer and [height
functions](heights.h).

We first overload the default function used to compute the normal,
defined in [fractions.h](). */
coord interface_normal (Point point, scalar c);
#define interface_normal(point, c) interface_normal (point, c)


#include "centered-weugene.h"
#include "rheology_model.h"
#include "tension.h"
#include "output_vtu_foreach.h"

/**
We will compute the normal using height-functions instead. If this is
not possible (typically at low resolutions) we revert back to
the Mixed-Youngs-Centered approximation. */

coord interface_normal (Point point, scalar c)
{
    coord n;
    if (!c.height.x.i || (n = height_normal (point, c, c.height)).x == nodata)
        n = mycs (point, c);
    return n;
}

//#define snapshot_i 5000
//#define dt_vtk 0.1
int snapshot_i = 500;
double dt_vtk = 0.05;
double Uin, Tin, Tcyl;
double RhoR, RhoRS, MuR, MuRS, CpR, CpRS, KappaR, KappaRS;
double Rho1, Rho2, Rho3;
double Mu0, Mu1, Mu2, Mu3;
double Kappa1, Kappa2, Kappa3;
double CP1, CP2, CP3;
double Ggrav, Ggrav_ndim;
double Sigma, sigma_ndim;
double cyl_diam, domain_size, dist_x, dist_y, front_x, Rbmin, Rbmax;
double ratio_Rbmin, ratio_Rbmax;
double ratio_dist_x, ratio_dist_y;
double ratio_front_x;
int Ncx, Ncy; //number of cylinders along Ox, Oy
int Nb; //number of bubbles
double Ca; // Ca = Mu*Ud/sigma
double Re; //Reynolds
double Fr; //Froude number Fr = sqrt(u^2/(g*cyl_diam))
double G;
double Umean;
double x_init = 2;
int maxlevel = 8;
int minlevel = 5;
int LEVEL = 9;
int adapt_method = 1; // 0 - traditional, 1 - using limitation, 2 - using array for maxlevel
double feps = 1e-10, fseps = 1e-10, ueps = 1e-3, Teps = 3e-2, aeps = 3e-2;
double TOLERANCE_P = 1e-6, TOLERANCE_V = 1e-7, TOLERANCE_T = 1e-6;

int main(int argc, char * argv[]) {
    TOLERANCE = 1e-7;
    NITERMIN = 1;
    NITERMAX = 100;
    CFL = 0.5;
    CFL_SIGMA = 0.7;
    CFL_ARR = 0.5;
    DT = 2e-5;
//	relative_residual_poisson = true;
// Physical parameters
	Uin = 1e-2, Tin = 300, Tcyl = 333.15;
//	Gorthala, R., Roux, J. A., & Vaughan, J. G. (1994). Resin Flow, Cure and Heat Transfer Analysis for Pultrusion Process. Journal of Composite Materials, 28(6), 486–506. doi:10.1177/002199839402800601  sci-hub.se/10.1177/002199839402800601
//EPON 9310
	Htr = 355000; //J/kg
	Arrhenius_const = 80600;//1/s
	Ea_by_R = 64000/8.314;// Kelvin
    n_degree = 1.2;

	Eeta_by_Rg = 3.76e+4/8.314;// Kelvin Epon 9310,Safonov page 21
	chi = 20;

	Rho1 = 1200, Rho2 = 1.092, Rho3 = 1790;//air at 23 C Graphite
//	Mu1 = 3.85e-7*exp(Eeta_by_Rg/Tin), Mu2 = 1e-4, Mu3 = Rho3*Mu1/Rho1/2.0;
    Mu0 = 3.85e-7, Mu1 = Mu0*exp(Eeta_by_Rg/Tin), Mu2 = 1.963e-5, Mu3 = Mu1;//air at 50C //Mu2 = 1.963e-5

	CP1 = 1255, CP2 = 1006, CP3 = 712;//J/(kg*K)
	Kappa1 = 0.2, Kappa2 = 0.02535, Kappa3 = 8.70;//W/(m*K)
	Sigma = 0.040;// N/m  0.040;
	Ggrav = 0; // m/s^2

	Nb = 10; Ncx = 5; Ncy = 7;
	ratio_Rbmin = 1./6.; ratio_Rbmax = 3./4.;
	ratio_dist_x = 2; ratio_dist_y = 2;
	ratio_front_x = -4;
	if (argc > 1)
		Tcyl = atof(argv[1]);
	if (argc > 2)
        maxlevel = atoi(argv[2]);
    //stokes = true;
	if (argc > 3)
        iter_fp = atoi(argv[3]);
	if (argc > 4)
        ratio_Rbmin = atof(argv[4]);
	if (argc > 5)
        ratio_Rbmax = atof(argv[5]);
	if (argc > 6)
        ratio_dist_x = atof(argv[6]);
	if (argc > 7)
        ratio_dist_y = atof(argv[7]);
	if (argc > 8 )
		ratio_front_x = atof(argv[8]);
	if (argc > 9)
        Nb = atoi(argv[9]);
	if (argc > 10)
        Ncx = atoi(argv[10]);
	if (argc > 11)
        Ncy = atoi(argv[11]);
	if (argc > 12)
        TOLERANCE_P = atof(argv[12]);
    if (argc > 13)
        TOLERANCE_V = atof(argv[13]);
    if (argc > 14)
        TOLERANCE_T = atof(argv[14]);
	if (argc > 15)
		Htr = atof(argv[15]);
	if (argc > 16)
		Arrhenius_const = atof(argv[16]);
	if (argc > 17)
		Ea_by_R = atof(argv[17]);
	cyl_diam = 15e-6; 
	dist_x = ratio_dist_x*cyl_diam, dist_y = ratio_dist_y*cyl_diam; //m 5-25 microns
	Rbmin = ratio_Rbmin*cyl_diam, Rbmax = ratio_Rbmax*cyl_diam;
	domain_size = dist_y*max(max(Ncx,Ncy),1);
	fprintf(ferr,
                 "Props(SI): Mu0=%g, Mu1=%g, Mu2=%g, Mu3=%g, Rho1=%g, Rho2=%g,  Rho3=%g,\n"
                 "           Kappa1=%g, Kappa2=%g, Kappa3=%g, CP1=%g, CP2=%g, CP3=%g,\n"
				 "           Sigma=%g, Uin=%g, time*=%g, Tin=%g, Tcyl=%g\n"
                 "           Htr=%g, Arrenius=%g, Ea_by_R=%g, n_deg=%g, m_deg=%g\n"
				 "           Eeta_by_Rg=%g, chi=%g\n"
                 "Apparatus: cyl_diam=%g,  domainSize=%g, Ncx=%d, Ncy=%d, Nb=%d\n",
                 Mu0, Mu1, Mu2, Mu3, Rho1, Rho2, Rho3,
                 Kappa1, Kappa2, Kappa3, CP1, CP2, CP3,
                 Sigma, Uin, cyl_diam/Uin, Tin, Tcyl,
                 Htr, Arrhenius_const, Ea_by_R, n_degree, m_degree,
                 Eeta_by_Rg, chi,
                 cyl_diam, domain_size, Ncx, Ncy, Nb);
// Dimensionless numbers
	Re = Uin*cyl_diam*Rho1/Mu1;
	Ca = Mu1*Uin/(Sigma + SEPS);
	Fr = sqrt(sq(Uin)/(Ggrav*cyl_diam + SEPS));
// Dimensionless parameters are chosen cyl_diam, rho1, Cp1, Tin, Uin
    L0 = domain_size/cyl_diam;
	front_x = ratio_front_x;
	Rbmin /= cyl_diam;
	Rbmax /= cyl_diam;
    dist_x /= cyl_diam;
	dist_y /= cyl_diam;
	RhoR = Rho2/Rho1, RhoRS = Rho3/Rho1;
	MuR = Mu2/Mu1, MuRS = Mu3/Mu1;
	CpR = CP2/CP1, CpRS = CP3/CP1;
	KappaR = Kappa2/Kappa1, KappaRS = Kappa3/Kappa1;
    rho1 = 1; rho2 = RhoR; rho3 = RhoRS;
    mu0 = (1./Re)*(Mu0/Mu1); mu1 = (1./Re); mu2 = mu1*MuR; mu3 = mu1*MuRS;
	Cp1 = 1; Cp2 = Cp1*CpR; Cp3 = Cp1*CpRS;//J/(kg*K)
	kappa1 = Kappa1/(Rho1*CP1*cyl_diam*Uin + SEPS), kappa2 = kappa1*KappaR, kappa3 = kappa1*KappaRS;//W/(m*K)
	Htr /= CP1*Tin;
	Arrhenius_const *= cyl_diam/(Uin + SEPS);
	Ea_by_R /= Tin;
	Eeta_by_Rg /= Tin;
	sigma_ndim = 1./(Re*Ca + SEPS);
    f.sigma = (fabs(sigma_ndim) < 1e+10) ? sigma_ndim : 0;
	Ggrav_ndim = 1./sq(Fr);
	Uin = 1;
	Tcyl /= Tin;
	Tin  /= Tin;
    const scalar temp_cyl[] = Tcyl;
    T_target = temp_cyl;
	cyl_diam = 1;
    origin (-L0/2, -L0/2.);
    N = 1 << LEVEL;
    periodic(top);
	fprintf(ferr,"Dim-less vars: mu0=%g, mu1=%g, mu2=%g, mu3=%g, rho1=%g, rho2=%g, rho3=%g,\n"
				 "               kappa1=%g, kappa2=%g, kappa3=%g, Cp1=%g, Cp2=%g, Cp3=%g,\n"
				 "               sigma=%g,  Uin=%g, Tin=%g, Tcyl=%g,\n"
				 "               Htr=%g, Arrhenius_const=%g, Ea_by_R=%g, Eeta_by_Rg=%g,\n"
				 "               L0=%g, cyl_diam=%g, dist_x=%g, dist_y=%g, front_x=%g, Rbmin=%g, Rbmax=%g,\n"
				 "               Ggrav_ndim=%g Uin=%g\n"
                 "Dim-less nums: Re=%g,  Ca=%g, Fr=%g\n"
                 "Solver:        DTmax=%g, CFL=%g, CFL_SIGMA=%g, CFL_ARR=%g, NITERMIN=%d,  NITERMAX=%d,\n"
                 "               TOLERANCE_P=%g, TOLERANCE_V=%g, TOLERANCE_T=%g\n"
                 "ADAPT:         minlevel=%d,  maxlevel=%d, feps=%g, fseps=%g, ueps=%g, Teps=%g, aeps=%g\n"
                 "OUTPUT:        dt_vtk=%g,    number of procs=%d\n",
				mu0, mu1, mu2, mu3, rho1, rho2, rho3,
				kappa1, kappa2, kappa3, Cp1, Cp2, Cp3,
				sigma_ndim, Uin, Tin, Tcyl,
				Htr, Arrhenius_const, Ea_by_R, Eeta_by_Rg,
				L0, cyl_diam, dist_x, dist_y, front_x, Rbmin, Rbmax,
				Ggrav_ndim, Uin,
                Re, Ca, Fr,
                DT, CFL, CFL_SIGMA, CFL_ARR, NITERMIN, NITERMAX,
                TOLERANCE_P, TOLERANCE_V, TOLERANCE_T,
                minlevel, maxlevel, feps, fseps, ueps, Teps, aeps,
                dt_vtk, npe());
    /**
    We must associate the height function field with the VOF tracer, so
    that it is used by the relevant functions (curvature calculation in
    particular). */
    f.height = h;
    fs.height = hs;
    a = av;
    fs.refine = fs.prolongation = fraction_refine;
    run();
}
//#define T_BC (0.5*(TMAX + TMIN) + 0.5*(TMAX - TMIN)*tanh((x)/(L0/10)))
//#define u_BC (Uin*(sq(0.5*L0) - sq(y)))
#define u_BC (1)
u.n[left]  = dirichlet(u_BC);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(1);
T[left]    = dirichlet(Tin);
fs[left]   = dirichlet(0);
alpha_doc[left] = dirichlet(0);//inflow is fresh resin

u.n[right] = neumann(0);
p[right]   = dirichlet(0);
pf[right]  = dirichlet(0);
T[right]    = neumann(0);
fs[right]   = neumann(0);
alpha_doc[right] = neumann(0);


#define RandMinMax(limMin, limMax) (0.5*((limMin)+(limMax)) + 0.5*noise()*((limMax) - (limMin)))
// ONLY for coord objects
#if dimension>2
  #define mynorm(v) (sqrt(sq(v.x) + sq(v.y) + sq(v.z)))
#else
  #define mynorm(v) (sqrt(sq(v.x) + sq(v.y)))
#endif

double sphere(double x, double y, double z, coord center, double radius) {
    return ( sq(x - center.x) + sq (y - center.y) + sq (z - center.z)
             - sq (radius));
}

double bubbles (double x, double y, double z)
{
    coord mypoint = {x,y,z};
    coord pnt_dist;
	double limMin = -L0/2. + Rbmax, limMax = front_x - Rbmax;
	double *R = malloc(Nb * sizeof(double));
    if (R == NULL) {
      fprintf(stderr, "malloc failed with R\n");
      return -1;
    }
	coord *centers = malloc(Nb * sizeof(coord));
	if (centers == NULL) {
      fprintf(stderr, "malloc failed with centers\n");
      return -1;
    }
// generating of Radii and centers of bubbles which are not overlapping
    srand (0);
    int iter = 0, i = 0;
	//fprintf(ferr, "bubble in\n");
    while(i < Nb){
        R[i] = RandMinMax(Rbmin, Rbmax);
        centers[i].x = RandMinMax(limMin, limMax);
        centers[i].y = RandMinMax(-L0/2.0 + Rbmax, L0/2.0 - Rbmax);
#if dimension>2
        centers[i].z = RandMinMax(-L0/2.0 + Rbmax, L0/2.0 - Rbmin);
#endif	
        for (int j = 0; j < i; j++) {
            foreach_dimension() pnt_dist.x = centers[i].x - centers[j].x;
            if ( mynorm(pnt_dist) < 1.3*(R[i] + R[j]) ) {
                i--;
                break; // toss again, because of overlapping
            };
        }
        i++;
        iter++;
        if (iter>100000*Nb) {fprintf(ferr, "too many fail attempts...\n"); exit(137);}
    }
// generating volume fraction f
    double phi = HUGE;
    for (int i = 0; i < Nb; i++) {
        //printf("i=%d x=%g y=%g R=%g\n", i, centers[i].x, centers[i].y, R[i] );
        foreach_dimension()
            pnt_dist.x = mypoint.x - centers[i].x;
        phi = min(phi, (mynorm(pnt_dist) - R[i]));
    }
	free(R);
	free(centers);
	//fprintf(ferr, "bubble end\n");
//    return min(phi, front_x - x);// with front
    return phi; // no front
}

double geometry(double x, double y, double z){
    coord mypoint = {x,y,z};
    coord pnt_dist;
	int k;
	double minv, maxv;
	double limMin = front_x + cyl_diam;
	double *R = malloc(Ncx*Ncy*sizeof(double));
    if (R == NULL) {
      fprintf(stderr, "malloc failed with R\n");
      return -1;
    }
	coord *centers = malloc(Ncx*Ncy*sizeof(coord));
	if (centers == NULL) {
      fprintf(stderr, "malloc failed with centers\n");
      return -1;
    }
//	fprintf(ferr, "fs in\n");
// generating of Radii and centers of bubbles which are not overlapping
    for (int i = 0; i < Ncx; i++){
		for(int j = 0; j < Ncy; j++){
			k = i*Ncy + j;
        	R[k] = RandMinMax(0.5*cyl_diam, 0.5*cyl_diam);
			minv = limMin + i*dist_x, maxv = limMin + i*dist_x;
        	centers[k].x = RandMinMax(minv, maxv);
			minv = -0.5*L0 + 0.5*dist_y + j*dist_y, maxv = -0.5*L0 + 0.5*dist_y + j*dist_y;
        	centers[k].y = RandMinMax(minv, maxv);
			//fprintf(ferr, "R=%g x=%g y=%g\n", R[k], centers[k].x, centers[k].y);
        }
    }
//	fprintf(ferr, "PHI\n");
// generating volume fraction f
    double phi = -HUGE;
	for (int i = 0; i < Ncx*Ncy; i++){
		foreach_dimension() pnt_dist.x = mypoint.x - centers[i].x;
       	phi = max(phi, (R[i] - mynorm(pnt_dist)));
//       	phi = -10;// delete
	}
	free(R);
	free(centers);
	//if (phi>0) fprintf(ferr, "phi=%g\n", phi);
//	fprintf(ferr, "fs out\n");
    return phi;

}

event init (t = 0) {
    if (!restore (file = "restart")) {
        int iter = 0;
        do {
            iter++;
			fprintf(ferr, "IN");
			fraction(fs, geometry(x,y,z));
			fraction(f, bubbles (x, y, z));
            filter_scalar(f, smoothed_f);
            filter_scalar(fs, smoothed_fs);
//            foreach() f[] = smoothed_f[];
//            foreach() fs[] = smoothed_fs[];
			fprintf(ferr, "ITER=%d\n", iter);
        }while ((iter <=3) || ((adapt_wavelet({smoothed_f, smoothed_fs}, (double []){feps, fseps}, maxlevel = maxlevel, minlevel=minlevel).nf != 0) && (iter <= 15)));
        fprintf(stderr, "init refinement iter=%d\n", iter);
        foreach() {
            T[] = Tin*(1 - fs[]) + fs[]*Tcyl;
            alpha_doc[] = 0;
            u.x[] = u_BC; //*(1-fs[]); // penalization will work
        }
        boundary ({T, alpha_doc, u});
        foreach_face(){
            double ff1 = (f[] + f[-1])/2.;
            double ff2 = (fs[] + fs[-1])/2.; //solid
            if (kappa1 || kappa2) {
                kappav.x[] = var_hom(ff1, ff2, kappa1, kappa2, kappa3);
            }
            if (mu1 || mu2) {
                double Tf = (T[] + T[-1])/2.;
                double alphaf = (alpha_doc[] + alpha_doc[-1])/2.;
                face vector muv = mu;
                muv.x[] = fm.x[] * mu(ff1, ff2, alphaf, Tf);
            }
        }
		boundary((scalar *){kappav, mu});
        heights (f, f.height);
        heights (fs, fs.height);
    }
	event("vtk_file");
}


event set_dtmax (i++) {
	RELATIVE_RES_TOLERANCE = 0.01;
    DT *= 1.05;
    DT = min(DT, 1e-2);
    fprintf(ferr, "set_dtmax: tnext= %g, t=%g, DT=%g, dt=%g\n", tnext, t, DT, dt);
}

event vof (i++) {
    if (f.height.x.i)
        heights (f, f.height);
    if (fs.height.x.i)
        heights (fs, fs.height);
}

/**
The gravity vector is aligned with the channel and viscosity is
unity. */

event acceleration (i++) {
	foreach_face(x)	av.x[] = Ggrav_ndim;
}

event advection_term(i++){
    TOLERANCE = TOLERANCE_P;
}

event viscous_term(i++){
    TOLERANCE = TOLERANCE_V;
    int m_bp = max(400 - 1*i, 1);
    double mindelta = L0/pow(2, maxlevel);
    double nu_bp = mu1/rho1;
    eta_s = sq(m_bp*mindelta)/nu_bp;
//    eta_s = 1e+6;
    fprintf(ferr, "m=%d, mindelta=%g, nu=%g, eta_s=%15.12g\n", m_bp, mindelta, nu_bp, eta_s);

}

event projection(i++){
    TOLERANCE = TOLERANCE_P;
}

event end_timestep(i++){
	TOLERANCE = TOLERANCE_T;
//	relative_residual_poisson = true;
}

//event vtk_file(i++){
//	relative_residual_poisson = false;
//}
event vtk_file (i+=500)
//event vtk_file (t += dt_vtk)
{
	char subname[80]; sprintf(subname, "heat_pol");
	scalar l[]; foreach() l[] = level;
	vector a_cell[]; foreach() foreach_dimension() a_cell.x[] = 0.5*(a.x[] + a.x[1]);
	vector mu_cell[]; foreach() foreach_dimension() mu_cell.x[] = 0.5*(mu.x[] + mu.x[1]);
//	vector uf_cell[]; foreach() foreach_dimension() uf_cell.x[] = 0.5*(uf.x[] + uf.x[1]);
//	vector kappa_cell[]; foreach() foreach_dimension() kappa_cell.x[] = 0.5*(kappa.x[] + kappa.x[1]);
    scalar curvature_by_sigma_cell[]; curvature (f, curvature_by_sigma_cell, f.sigma, add = false);
//    f.height = h;
//    fs.height = hs;
    if (f.height.x.i)
        heights (f, f.height);
    if (fs.height.x.i)
        heights (fs, fs.height);
    vector nnf[], nnfs[];
    foreach() {
//        coord nnfc = interface_normal (point, f);
//        coord nnfsc = interface_normal (point, fs);
        coord nnfc = height_normal (point, f, f.height);
        coord nnfsc = height_normal (point, fs, fs.height);
        foreach_dimension(){
            nnf.x[] = nnfc.x;
            nnfs.x[] = nnfsc.x;
        }
    }
    boundary((scalar *){nnf, nnfs});
#ifdef DEBUG_BRINKMAN_PENALIZATION
	output_vtu_MPI(subname, (iter_fp) ? t + dt : 0, (scalar *) {T, alpha_doc, p, fs, f, l,  rhov}, (vector *) {u, dbp, total_rhs, av, residual_of_u, conv_term, mu, kappa});
#else
	fprintf(ferr, "output_vtu_MPI");
//    output_vtu_MPI(subname, (iter_fp) ? t + dt : 0, (scalar *) {T, alpha_doc, p, fs, f}, (vector *) {u});
    output_vtu_MPI(subname, (iter_fp) ? t + dt : 0, list = (scalar *) {T, alpha_doc, p, fs, f, l, rhov, smoothed_f, smoothed_fs, curvature_by_sigma_cell},
                   vlist = (vector *) {u, a_cell, mu_cell, g, nnf, nnfs, av});
#endif
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

//#define ADAPT_SCALARS {f, fs, un, T, alpha_doc}
#define ADAPT_SCALARS {smoothed_f, smoothed_fs, un, T, alpha_doc}
#define ADAPT_EPS_SCALARS {feps, fseps, ueps, Teps, aeps}
//#define ADAPT_SCALARS {f, fs}
//#define ADAPT_EPS_SCALARS {feps, fseps}

scalar un[]; 
event adapt (i++){
	foreach() un[] = norm(u); boundary((scalar *){un});
	double eps_arr[] = ADAPT_EPS_SCALARS;
	MinMaxValues((scalar *) ADAPT_SCALARS, eps_arr);
	adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
	fraction(fs, geometry(x,y,z));
	//if (i > 300) stokes = true;
}

event stop(t = 100 * L0 / Uin);
