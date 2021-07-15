#define BRINKMAN_PENALIZATION 1
#define DEBUG_MINMAXVALUES
//#define DEBUG_BRINKMAN_PENALIZATION
#define DEBUG_MODE_POISSON
#define REACTION_MODEL REACTION_MODEL_NON_AUTOCATALYTIC
//#define DEBUG_OUTPUT_VTU_MPI
#define FILTERED
#define JACOBI 1
//#define RELATIVE_RES
//#define PRINT_ALL_VALUES
//#define STOKES
#define DIRICHLET_BC 0
scalar omega[];
scalar l2[], un[];
face vector fs_face[];
//face vector av[];


#include "centered-weugene.h"
#include "rheology_model.h"
#include "tension.h"
#include "output_vtu_foreach.h"
#include "tag.h"
#include <ctype.h>
char dump_name[30];

int snapshot_i = 5000;
double dt_vtk = 0.1;
double Uin, Tin, Tcyl;
double RhoR, RhoRS, MuR, MuRS, CpR, CpRS, KappaR, KappaRS;
double Rho1, Rho2, Rho3;
double Mu0, Mu1, Mu2, Mu3;
double Kappa1, Kappa2, Kappa3;
double CP1, CP2, CP3;
double Ggrav, Ggrav_ndim;
double Sigma, sigma_ndim;
double cyl_diam, domain_size, dist_x, dist_y, cyl_x, front_x, Rbmin, Rbmax;
double ratio_Rbmin, ratio_Rbmax;
double ratio_dist_x, ratio_dist_y, dev_r, develx, devely;
double ratio_front_x;
double shift_x, shift_y;
int non_saturated;
int Ncx, Ncy; //number of cylinders along Ox, Oy
int Nb; //number of bubbles
double Ca; // Ca = Mu*Ud/sigma
double Re; //Reynolds
double Fr; //Froude number Fr = sqrt(u^2/(g*cyl_diam))
double G;
double Umean;
double x_init = 2, Dx_min, dx_min;
int maxlevel = 9;
int minlevel = 5;
int LEVEL = 7;
double myt = 0;
int adapt_method = 1; // 0 - traditional, 1 - using limitation, 2 - using array for maxlevel
double feps = 1e-10, fseps = 1e-10, ueps = 1e-2, Teps = 3e-2, aeps = 3e-2;
double TOLERANCE_P = 1e-6, TOLERANCE_V = 1e-7, TOLERANCE_T = 1e-6;
double *R = NULL;
coord *centers = NULL;

//get time from dump name (dump-1.1234 => 1.1234)
double get_double(const char *str)
{
    /* First skip non-digit characters */
    /* Special case to handle negative numbers and the `+` sign */
    while (*str && !(isdigit(*str) || ((*str == '-' || *str == '+') && isdigit(*(str + 1)))))
        str++;
    /* The parse to a double */
    return fabs(strtod(str, NULL));
}

int main(int argc, char * argv[]) {
    TOLERANCE = 1e-7;
    NITERMIN = 1;
    NITERMAX = 100;
    CFL = 0.5;
    CFL_SIGMA = 0.3;
    CFL_ARR = 0.5;
    DT = 1e-5;
    m_bp = 1; // Brinkman layer resolution
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
    CP1 = 1255, CP2 = 1006, CP3 = 712;//J/(kg*K)
    Kappa1 = 0.2, Kappa2 = 0.02535, Kappa3 = 8.70;//W/(m*K)
    Sigma = 0.040*1e-1;// N/m  0.040;
    Ggrav = 0; // m/s^2

    Nb = 10; Ncx = 5; Ncy = 8;
    ratio_Rbmin = 1./6.; ratio_Rbmax = 3./4.;
    ratio_dist_x = 2; ratio_dist_y = 2;
    ratio_front_x = -6;
    cyl_x = -5;
    shift_x = 0;
    shift_y = 0;
    dev_r = 0.4, develx = 0.4, devely = 0.4;
    non_saturated = 1;

    // set which dump files will be converted: each $(i_take)th
    // by default each dump will be converted
    if (argc > 1)
        strcpy(dump_name, argv[1]);
    else
        return 1;
    if (argc > 2)
        iter_fp = atoi (argv[2]);
    if (argc > 3)
        maxlevel = atoi (argv[3]);
    //	Mu1 = 3.85e-7*exp(Eeta_by_Rg/Tin), Mu2 = 1e-4, Mu3 = Rho3*Mu1/Rho1/2.0;
    Mu0 = 3.85e-7, Mu1 = Mu0*exp(Eeta_by_Rg/Tin), Mu2 = 1.963e-5, Mu3 = Mu1*Rho3/Rho1;//air at 50C //Mu2 = 1.963e-5

    myt = get_double(dump_name);
    fprintf(ferr, "dump_name=%s iter_fp=%d maxlevel=%d\n", dump_name, iter_fp, maxlevel);

	cyl_diam = 15e-6;
	dist_x = ratio_dist_x*cyl_diam, dist_y = ratio_dist_y*cyl_diam; //m 5-25 microns
	Rbmin = ratio_Rbmin*cyl_diam, Rbmax = ratio_Rbmax*cyl_diam;
	domain_size = dist_y*max(max(Ncx,Ncy),1);
    Dx_min = domain_size/pow(2.0, maxlevel);
//    Mu3 = Mu1*(Rho3/Rho1);
	fprintf(ferr,
                 "Props(SI): Mu0=%g, Mu1=%g, Mu2=%g, Mu3=%g, Rho1=%g, Rho2=%g,  Rho3=%g,\n"
                 "           Kappa1=%g, Kappa2=%g, Kappa3=%g, CP1=%g, CP2=%g, CP3=%g,\n"
				 "           Sigma=%g, Uin=%g, time*=%g, Tin=%g, Tcyl=%g\n"
                 "           Htr=%g, Arrenius=%g, Ea_by_R=%g, n_deg=%g, m_deg=%g\n"
				 "           Eeta_by_Rg=%g, chi=%g\n"
                 "Apparatus: cyl_diam=%g,  domainSize=%g, Dx_min=%g, Ncx=%d, Ncy=%d, Nb=%d\n",
                 Mu0, Mu1, Mu2, Mu3, Rho1, Rho2, Rho3,
                 Kappa1, Kappa2, Kappa3, CP1, CP2, CP3,
                 Sigma, Uin, cyl_diam/Uin, Tin, Tcyl,
                 Htr, Arrhenius_const, Ea_by_R, n_degree, m_degree,
                 Eeta_by_Rg, chi,
                 cyl_diam, domain_size, Dx_min, Ncx, Ncy, Nb);
	Ncy += (shift_y > 0); //for saturated flow we add artificial cylinder
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
    dx_min = L0/pow(2, maxlevel);
	RhoR = Rho2/Rho1, RhoRS = Rho3/Rho1;
	MuR = Mu2/Mu1, MuRS = Mu3/Mu1;
	CpR = CP2/CP1, CpRS = CP3/CP1;
	KappaR = Kappa2/Kappa1, KappaRS = Kappa3/Kappa1;
    rho1 = 1; rho2 = RhoR; rho3 = RhoRS;
    mu0 = (1./Re)*(Mu0/Mu1); mu1 = (1./Re); mu2 = mu1*MuR; mu3 = mu1*MuRS;
//    mu3 = mu1*(rho3/rho1)*sq(1.0/(m_bp*dx_min));
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
                 "               dev_r=%g develx=%g devely=%g ldomain=%g\n"
                 "               shift_x=%g shift_y=%g non_saturated=%d\n"
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
                dev_r, develx, devely, domain_size/L0,
                shift_x, shift_y, non_saturated,
                DT, CFL, CFL_SIGMA, CFL_ARR, NITERMIN, NITERMAX,
                TOLERANCE_P, TOLERANCE_V, TOLERANCE_T,
                minlevel, maxlevel, feps, fseps, ueps, Teps, aeps,
                dt_vtk, npe());
//    a = av;
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


event init (t = 0) {
    bool success = restore (file = dump_name);
    fprintf(ferr, "file has been read: L0=%g\n", L0);

    if (!success) {
        fprintf(ferr, "can't open the file %s. Missing this file, go to the next file\n", dump_name);
        return 0;
    }
    event("vtk_file");
    exit(777);
}

event vtk_file (i++)
{
    char subname[80]; sprintf(subname, "dump2pvd_heat");
    scalar l[]; foreach() l[] = level;
    scalar mu_cell[], dpdx[];
    foreach() {
        double mus = 0;
        foreach_dimension() {
            mus += mu.x[] + mu.x[1];
        }
        mu_cell[] = mus/(2.0*dimension);
        dpdx[] = (p[1] - p[-1])/(2*Delta);
    }
    boundary((scalar *){mu_cell});

    fprintf(ferr, "output_vtu_MPI is running...\n");
    output_vtu_MPI(subname, myt, list = (scalar *) {T, dpdx, alpha_doc, p, fs, f, l, rhov, mu_cell},
    vlist = (vector *) {u, a});//, mu_cell_minus, mu_cell_plus
    return 1;
}

event stop(t = 10*L0 / Uin);
