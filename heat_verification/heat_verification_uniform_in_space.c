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


#include "centered-weugene.h"
#include "rheology_model.h"
#include "output_vtu_foreach.h"

//#define snapshot_i 5000
//#define dt_vtk 0.1
int snapshot_i = 1000;
double dt_vtk = 0.1;
double Uin, Tin, Tcyl;
double RhoR, RhoRS, MuR, MuRS, CpR, CpRS, KappaR, KappaRS;
double Rho1, Rho2, Rho3;
double Mu0, Mu1, Mu2, Mu3;
double Kappa1, Kappa2, Kappa3;
double CP1, CP2, CP3;
double domain_size;
double Re; //Reynolds
double Umean;
int maxlevel = 7;
int minlevel = 5;
int LEVEL = 7;
int adapt_method = 1; // 0 - traditional, 1 - using limitation, 2 - using array for maxlevel
double feps = 1e-10, fseps = 1e-10, ueps = 1e-3, Teps = 3e-2, aeps = 3e-2;
double TOLERANCE_P = 1e-6, TOLERANCE_V = 1e-7, TOLERANCE_T = 1e-6;

int main(int argc, char * argv[]) {
    TOLERANCE = 1e-7;
    NITERMIN = 1;
    NITERMAX = 100;
    CFL = 0.1;
    CFL_ARR = 0.1;
    DT = 2e-5;
// Physical parameters
	Uin = 1e-2, Tcyl = Tin = 400, Tcyl = 333.15;
//	Gorthala, R., Roux, J. A., & Vaughan, J. G. (1994). Resin Flow, Cure and Heat Transfer Analysis for Pultrusion Process. Journal of Composite Materials, 28(6), 486–506. doi:10.1177/002199839402800601  sci-hub.se/10.1177/002199839402800601
//EPON 9310
	Htr = 355000; //J/kg
	Arrhenius_const = 80600;//1/s
	Ea_by_R = 64000/8.314;// Kelvin

	Eeta_by_Rg = 3.76e+4/8.314;// Kelvin Epon 9310, Safonov page 21
	chi = 20;
	n_degree = 1.2;

	Rho1 = 1200, Rho2 = 10.092, Rho3 = 2560;//air at 50C
//	Mu1 = 3.85e-7*exp(Eeta_by_Rg/Tin), Mu2 = 1e-4, Mu3 = Rho3*Mu1/Rho1/2.0;
    Mu0 = 3.85e-7, Mu1 = Mu0*exp(Eeta_by_Rg/Tin), Mu2 = 1.963e-5, Mu3 = Mu1;//air at 50C //Mu2 = 1.963e-5

	CP1 = 1255, CP2 = 1006, CP3 = 670;//J/(kg*K)
	Kappa1 = 0.2, Kappa2 = 0.02535, Kappa3 = 1.04;//W/(m*K)

	if (argc > 1)
        maxlevel = atoi(argv[1]);
    //stokes = true;
	if (argc > 2)
        iter_fp = atoi(argv[2]);
	if (argc > 3)
        TOLERANCE_P = atof(argv[3]);
    if (argc > 4)
        TOLERANCE_V = atof(argv[4]);
    if (argc > 5)
        TOLERANCE_T = atof(argv[5]);
	if (argc > 6)
		Htr = atof(argv[6]);
	if (argc > 7)
		Arrhenius_const = atof(argv[7]);
	if (argc > 8)
		Ea_by_R = atof(argv[8]);
	domain_size = 1e-2;
	fprintf(ferr,
                 "Props(SI): Mu0=%g, Mu1=%g, Mu2=%g, Mu3=%g, Rho1=%g, Rho2=%g,  Rho3=%g,\n"
                 "           Kappa1=%g, Kappa2=%g, Kappa3=%g, CP1=%g, CP2=%g, CP3=%g,\n"
				 "           Uin=%g, time*=%g, Tin=%g, Tcyl=%g,\n"
                 "           Htr=%g, Arrenius=%g, Ea_by_R=%g, n_deg=%g, m_deg=%g,\n"
				 "           Eeta_by_Rg=%g, chi=%g\n"
                 "Apparatus: domainSize=%g\n",
                 Mu0, Mu1, Mu2, Mu3, Rho1, Rho2, Rho3,
                 Kappa1, Kappa2, Kappa3, CP1, CP2, CP3,
                 Uin, domain_size/Uin, Tin, Tcyl,
                 Htr, Arrhenius_const, Ea_by_R, n_degree, m_degree,
                 Eeta_by_Rg, chi,
                 domain_size);
// Dimensionless numbers
	Re = Uin*domain_size*Rho1/Mu1;
// Dimensionless parameters are chosen cyl_diam, rho1, Cp1, Tin, Uin
    L0 = 1;
	RhoR = Rho2/Rho1, RhoRS = Rho3/Rho1;
	MuR = Mu2/Mu1, MuRS = Mu3/Mu1;
	CpR = CP2/CP1, CpRS = CP3/CP1;
	KappaR = Kappa2/Kappa1, KappaRS = Kappa3/Kappa1;
    rho1 = 1; rho2 = RhoR; rho3 = RhoRS;
    mu0 = (1./Re)*(Mu0/Mu1); mu1 = (1./Re); mu2 = mu1*MuR; mu3 = mu1*MuRS;
	Cp1 = 1; Cp2 = Cp1*CpR; Cp3 = Cp1*CpRS;//J/(kg*K)
	kappa1 = Kappa1/(Rho1*CP1*domain_size*Uin + SEPS), kappa2 = kappa1*KappaR, kappa3 = kappa1*KappaRS;//W/(m*K)
	Htr /= CP1*Tin;
	Arrhenius_const *= domain_size/(Uin + SEPS);
	Ea_by_R /= Tin;
	Eeta_by_Rg /= Tin;
	Uin = 1;
	Tcyl /= Tin;
	Tin  /= Tin;
    domain_size = 1;
    origin (-L0/2, -L0/2.);
    N = 1 << LEVEL;
    periodic(right);
    periodic(top);
    stokes = true;
    stokes_heat = true;
	fprintf(ferr,"Dim-less vars: mu0=%g, mu1=%g, mu2=%g, mu3=%g, rho1=%g, rho2=%g, rho3=%g,\n"
				 "               kappa1=%g, kappa2=%g, kappa3=%g, Cp1=%g, Cp2=%g, Cp3=%g,\n"
				 "               Uin=%g, Tin=%g,\n"
				 "               Htr=%g, Arrhenius_const=%g, Ea_by_R=%g, Eeta_by_Rg=%g,\n"
				 "               L0=%g, Uin=%g\n"
                 "Dim-less nums: Re=%g\n"
                 "Solver:        DTmax=%g, CFL=%g, CFL_ARR=%g, NITERMIN=%d,  NITERMAX=%d,\n"
                 "               TOLERANCE_P=%g, TOLERANCE_V=%g, TOLERANCE_T=%g stokes=%d stokes_heat=%d\n"
                 "ADAPT:         minlevel=%d,  maxlevel=%d, feps=%g, fseps=%g, ueps=%g, Teps=%g, aeps=%g\n"
                 "OUTPUT:        dt_vtk=%g,    number of procs=%d\n",
				mu0, mu1, mu2, mu3, rho1, rho2, rho3,
				kappa1, kappa2, kappa3, Cp1, Cp2, Cp3,
				Uin, Tin,
				Htr, Arrhenius_const, Ea_by_R, Eeta_by_Rg,
				L0, Uin,
                Re,
                DT, CFL, CFL_ARR, NITERMIN, NITERMAX,
                TOLERANCE_P, TOLERANCE_V, TOLERANCE_T, stokes, stokes_heat,
                minlevel, maxlevel, feps, fseps, ueps, Teps, aeps,
                dt_vtk, npe());
    fs.refine = fs.prolongation = fraction_refine;
    run();
}
//#define T_BC (0.5*(TMAX + TMIN) + 0.5*(TMAX - TMIN)*tanh((x)/(L0/10)))
//#define u_BC (Uin*(sq(0.5*L0) - sq(y)))
#define u_BC (0)

// ONLY for coord objects
#if dimension>2
  #define mynorm(v) (sqrt(sq(v.x) + sq(v.y) + sq(v.z)))
#else
  #define mynorm(v) (sqrt(sq(v.x) + sq(v.y)))
#endif


event init (t = 0) {
    if (!restore (file = "restart")) {

        foreach() {
            fs[] = 0;
            f[] = 1;
            T[] = Tin*(1 - fs[]) + fs[]*Tcyl;
            alpha_doc[] = 0;
            u.x[] = 0; //*(1-fs[]); // penalization will work
        }
        boundary ({f, fs, T, alpha_doc, u});
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
    }
	event("vtk_file");
}


event set_dtmax (i++) {
	RELATIVE_RES_TOLERANCE = 0.01;
    DT *= 1.05;
    DT = min(DT, 1e-3);
    fprintf(ferr, "set_dtmax: tnext= %g, t=%g, DT=%g, dt=%g\n", tnext, t, DT, dt);
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


//event vtk_file (i+=500)
event vtk_file (t += dt_vtk)
{
	char subname[80]; sprintf(subname, "heat_pol");
	scalar l[]; foreach() l[] = level;
#ifdef DEBUG_BRINKMAN_PENALIZATION
	output_vtu_MPI(subname, (iter_fp) ? t + dt : 0, (scalar *) {T, alpha_doc, p, l,  rhov}, (vector *) {u, dbp, total_rhs, residual_of_u, mu, kappa});
#else
	fprintf(ferr, "output_vtu_MPI");
//    output_vtu_MPI(subname, (iter_fp) ? t + dt : 0, (scalar *) {T, alpha_doc, p, fs, f}, (vector *) {u});
    output_vtu_MPI(subname, (iter_fp) ? t + dt : 0, list = (scalar *) {T, alpha_doc, p, l, fs, f},
                   vlist = (vector *) {u});
#endif
}


#if DUMP
event snapshot (i += snapshot_i)
//event snapshot (t += 1e-1)
{
    char name[80];
    sprintf(name, "dump-%04g",t);
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

//scalar un[];
//event adapt (i++){
//	foreach() un[] = norm(u); boundary((scalar *){un});
//	double eps_arr[] = ADAPT_EPS_SCALARS;
//	MinMaxValues((scalar *) ADAPT_SCALARS, eps_arr);
//	adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
//	fraction(fs, geometry(x,y,z));
//	//if (i > 300) stokes = true;
//}

event stop(t = 100*L0/Uin);