#define BRINKMAN_PENALIZATION 1
//#define DEBUG_MINMAXVALUES
//#define DEBUG_BRINKMAN_PENALIZATION
//#define DEBUG_MODE_TENSION
//#define DEBUG_MODE_POISSON
//#define DEBUG_HEAT
#define REACTION_MODEL REACTION_MODEL_NON_AUTOCATALYTIC
//#define DEBUG_OUTPUT_VTU_MPI
#define FILTERED
#define JACOBI 1
#define CORRECT_UF_FLUXES 1
//#define CURVATURE_CORR // see tension.h
#define RELATIVE_RES
#define EPS_MAXA 2 // adapt based on Max-Min value
//#define PRINT_ALL_VALUES
#ifdef DEBUG_MODE_TENSION
scalar f_hat[];
#endif
//#define STOKES
#define T_DIRICHLET_BC 1
#define STICKY_SOLID 1
scalar which_meth[];
scalar un[], mu_cell[];
face vector fs_face[];
static coord vel_s = {0, 0, 0};
//#include "curvature_partstr.h"
#include "centered-weugene.h"
#include "rheology_model.h"
#include "tension.h"
#include "utils-weugene.h"
#include "output_vtu_foreach.h"
#include "tag.h"

int snapshot_i = 1000;
double dt_vtk = 0.1;
double Uin, Tin, Tcyl, Tam;
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
double maxDT;
int adapt_method = 1; // 0 - traditional, 1 - using limitation, 2 - using array for maxlevel
double feps = 1e-10, fseps = 1e-10, ueps = 1e-2, rhoeps = 1e-10, Teps = 3e-2, aeps = 3e-2, mueps=1e-2;
double TOLERANCE_P = 1e-7, TOLERANCE_V = 1e-8, TOLERANCE_T = 1e-7;
double rslice = 1.2;
char subname[150], logname[80];
int main(int argc, char * argv[]) {
    fprintf(ferr,"correction of uf\n");
    fprintf(ferr, "./a.out Tcyl, Tin, maxlevel, iter_fp, ratio_Rbmin, ratio_Rbmax, ratio_dist_x, ratio_dist_y, ratio_front_x, "
                  "cyl_x, Nb, Ncx, Ncy, TOLERANCE_P, TOLERANCE_V, TOLERANCE_T, Htr, Arrhenius_const, Ea_by_R, dev_r, develx, devely, "
                  "shift_x, shift_y, non_saturated\n");
    TOLERANCE = 1e-7;
    NITERMIN = 1;
    NITERMAX = 100;
    CFL = 0.4;
    ignore_tension_CFL = false;
    CFL_ARR = 0.5;
    DT = 1e-6;
    N_smooth = 1; //three-phase-rheology.h
#ifdef STOKES
    stokes = false;
    stokes_heat = false;
#endif
    m_bp = 1; // Brinkman layer resolution, then it reduces to 0.1
    fprintf(ferr, "Stokes flow:stokes=%d, stokes_heat=%d, uf.x correction. Rho2 is 10x and Mu2 is 100x higher\n", stokes, stokes_heat);
#if T_DIRICHLET_BC != 0
    fprintf(ferr, "T_DIRICHLET_BC: 1\n");
#endif
    sprintf(subname, "saturated_1");
    sprintf(logname, "log_saturated_1");
// Physical parameters
	Uin = 1e-2, Tin = 300, Tcyl = 400.15, Tam = 300;
//	Gorthala, R., Roux, J. A., & Vaughan, J. G. (1994). Resin Flow, Cure and Heat Transfer Analysis for Pultrusion Process. Journal of Composite Materials, 28(6), 486–506. doi:10.1177/002199839402800601  sci-hub.se/10.1177/002199839402800601
//EPON 9310
	Htr = 355000; //J/kg
	Arrhenius_const = 80600;//1/s //originally in Gothala: 80600
	Ea_by_R = 12000/8.314;// Kelvin //origianally  64000/8.314;
    n_degree = 1.2;

	Eeta_by_Rg = 3.76e+4/8.314;// Kelvin Epon 9310,Safonov page 21
	chi = 20;
    Rho1 = 1200, Rho2 = 1.092, Rho3 = 1790;//air at 23 C Graphite
    CP1 = 1255, CP2 = 1006, CP3 = 712;//J/(kg*K)
    Kappa1 = 0.2, Kappa2 = 0.02535, Kappa3 = 8.70;//W/(m*K)
    Sigma = 0.040;// N/m  0.040;
    Ggrav = 0; // m/s^2

    Nb = 10; Ncx = 5; Ncy = 8;
    ratio_Rbmin = 1./6.; ratio_Rbmax = 1.2;
    ratio_dist_x = 2; ratio_dist_y = 2;
    ratio_front_x = -6;
    cyl_x = -5;
    shift_x = 0;
    shift_y = 0;
    dev_r = 0.4, develx = 0.4, devely = 0.4;
    non_saturated = 1;


    if (argc > 1)
        Tcyl = atof(argv[1]);
    if (argc > 2)
        Tin = atof(argv[2]);
//	Mu1 = 3.85e-7*exp(Eeta_by_Rg/Tin), Mu2 = 1e-4;
    Mu0 = 3.85e-7, Mu1 = Mu0*exp(Eeta_by_Rg/Tin), Mu2 = 1.963e-5, Mu3 = Mu1*Rho3/Rho1;//air at 50C //Mu2 = 1.963e-5
	if (argc > 3)
        maxlevel = atoi(argv[3]);
    if (argc > 4)
        iter_fp = atoi(argv[4]);
	if (argc > 5)
        ratio_Rbmin = atof(argv[5]);
	if (argc > 6)
        ratio_Rbmax = atof(argv[6]);
	if (argc > 7)
        ratio_dist_x = atof(argv[7]);
	if (argc > 8)
        ratio_dist_y = atof(argv[8]);
	if (argc > 9)
		ratio_front_x = atof(argv[9]);
    if (argc > 10 )
        cyl_x = atof(argv[10]);
	if (argc > 11)
        Nb = atoi(argv[11]);
	if (argc > 12)
        Ncx = atoi(argv[12]);
	if (argc > 13)
        Ncy = atoi(argv[13]);
	if (argc > 14)
        TOLERANCE_P = atof(argv[14]);
    if (argc > 15)
        TOLERANCE_V = atof(argv[15]);
    if (argc > 16)
        TOLERANCE_T = atof(argv[16]);
	if (argc > 17)
		Htr = atof(argv[17]);
	if (argc > 18)
		Arrhenius_const = atof(argv[18]);
	if (argc > 19)
		Ea_by_R = atof(argv[19]);
    if (argc > 20)
        dev_r = atof(argv[20]);
    if (argc > 21)
        develx = atof(argv[21]);
    if (argc > 22)
        devely = atof(argv[22]);
    if (argc > 23)
        shift_x = atof(argv[23]);
    if (argc > 24)
        shift_y = atof(argv[24]);
    if (argc > 25)
        non_saturated = atoi(argv[25]);
    if (argc > 26) {
        strcpy(subname, argv[26]);
        sprintf (logname, "log%s", subname);
    }
    maxDT = (non_saturated>0) ? 1e-3 : 1e-2;//5e-04;
	cyl_diam = 15e-6;
	dist_x = ratio_dist_x*cyl_diam, dist_y = ratio_dist_y*cyl_diam; //m 5-25 microns
	Rbmin = ratio_Rbmin*cyl_diam, Rbmax = ratio_Rbmax*cyl_diam;
	domain_size = dist_y*max(max(Ncx,Ncy),1);
    Dx_min = domain_size/pow(2.0, maxlevel);
	fprintf(ferr,
                 "Props(SI): Mu0=%g, Mu1=%g, Mu2=%g, Mu3=%g, Rho1=%g, Rho2=%g,  Rho3=%g,\n"
                 "           nu1=%g, nu2=%g, nu3=%g,\n"
                 "           Kappa1=%g, Kappa2=%g, Kappa3=%g, CP1=%g, CP2=%g, CP3=%g,\n"
				 "           Sigma=%g, Uin=%g, time*=%g, Tin=%g, Tcyl=%g, Tam=%g\n"
                 "           Htr=%g, Arrenius=%g, Ea_by_R=%g, n_deg=%g, m_deg=%g\n"
				 "           Eeta_by_Rg=%g, chi=%g\n"
                 "Geometry: cyl_diam=%g,  domainSize=%g, Dx_min=%g, Ncx=%d, Ncy=%d, Nb=%d\n",
                 Mu0, Mu1, Mu2, Mu3, Rho1, Rho2, Rho3,
                 Mu1/Rho1, Mu2/Rho2, Mu3/Rho3,
                 Kappa1, Kappa2, Kappa3, CP1, CP2, CP3,
                 Sigma, Uin, cyl_diam/Uin, Tin, Tcyl, Tam,
                 Htr, Arrhenius_const, Ea_by_R, n_degree, m_degree,
                 Eeta_by_Rg, chi,
                 cyl_diam, domain_size, Dx_min, Ncx, Ncy, Nb);
	Ncy += (shift_y > 0); //for saturated flow we add artificial cylinder
// Dimensionless numbers
	Re = Uin*cyl_diam*Rho1/Mu1;
	Ca = Mu1*Uin/(Sigma + SEPS);
	Fr = sqrt(sq(Uin)/(Ggrav*cyl_diam + SEPS));
// Dimensionless parameters are chosen cyl_diam, rho1, Cp1, Tin, Uin
    size(20);
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
    mu0 = (1./Re)*(Mu0/Mu1); mu1 = (1./Re); mu2 = mu1*MuR; mu3 = mu1*MuRS; mu_eff = 0;
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
    Tam /= Tin;
	Tcyl /= Tin;
	Tin  /= Tin;
    const scalar temp_cyl[] = Tcyl;
    T_target = temp_cyl;
    const vector U_sol[] = {vel_s.x, vel_s.y, vel_s.z};
    target_U = U_sol;
	cyl_diam = 1;
    init_grid(1 << LEVEL);
    origin (-L0/2, -L0/2.);
//    periodic(right);
//    periodic(top);
	fprintf(ferr,"Dim-less vars: mu0=%g, mu1=%g, mu2=%g, mu3=%g, rho1=%g, rho2=%g, rho3=%g,\n"
                 "               nu1=%g, nu2=%g, nu3=%g,\n"
				 "               kappa1=%g, kappa2=%g, kappa3=%g, Cp1=%g, Cp2=%g, Cp3=%g,\n"
				 "               sigma=%g,  Uin=%g, Tin=%g, Tcyl=%g, Tam=%g,\n"
				 "               Htr=%g, Arrhenius_const=%g, Ea_by_R=%g, Eeta_by_Rg=%g,\n"
				 "               L0=%g, cyl_diam=%g, dist_x=%g, dist_y=%g, front_x=%g, Rbmin=%g, Rbmax=%g,\n"
				 "               Ggrav_ndim=%g Uin=%g\n"
                 "Dim-less nums: Re=%g,  Ca=%g, Fr=%g\n"
                 "               dev_r=%g develx=%g devely=%g ldomain=%g\n"
                 "               shift_x=%g shift_y=%g non_saturated=%d\n"
                 "Solver:        DTmax=%g, CFL=%g, CFL_ARR=%g, NITERMIN=%d,  NITERMAX=%d,\n"
                 "               TOLERANCE_P=%g, TOLERANCE_V=%g, TOLERANCE_T=%g\n"
                 "ADAPT:         minlevel=%d,  maxlevel=%d, feps=%g, fseps=%g, ueps=%g, Teps=%g, aeps=%g\n"
                 "OUTPUT:        dt_vtk=%g,    number of procs=%d\n",
				mu0, mu1, mu2, mu3, rho1, rho2, rho3,
				mu1/rho1, mu2/rho2, mu3/rho3,
				kappa1, kappa2, kappa3, Cp1, Cp2, Cp3,
				sigma_ndim, Uin, Tin, Tcyl, Tam,
				Htr, Arrhenius_const, Ea_by_R, Eeta_by_Rg,
				L0, cyl_diam, dist_x, dist_y, front_x, Rbmin, Rbmax,
				Ggrav_ndim, Uin,
                Re, Ca, Fr,
                dev_r, develx, devely, domain_size/L0,
                shift_x, shift_y, non_saturated,
                DT, CFL, CFL_ARR, NITERMIN, NITERMAX,
                TOLERANCE_P, TOLERANCE_V, TOLERANCE_T,
                minlevel, maxlevel, feps, fseps, ueps, Teps, aeps,
                dt_vtk, npe());
//    a = av;
    fs.refine = fs.prolongation = fraction_refine;
    f.refine = f.prolongation = fraction_refine;
    boundary((scalar *){f, fs});
#ifdef _MPI
    int rank, psize, h_len;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    // get rank of this proces
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // get total process number
    MPI_Comm_size(MPI_COMM_WORLD, &psize);
    MPI_Get_processor_name(hostname, &h_len);
    printf("rank:%d size: %d at %s h_len %d\n", rank, psize, hostname, h_len);
#endif
#if CURV_PARTSTR==1
    DumpCsvInit();
    partstr_conf.nohf = false; // both GHF and particle curvature are possible
#endif
    run();
#if CURV_PARTSTR==1
    DumpCsvFin();
#endif
}
//#define T_BC (0.5*(TMAX + TMIN) + 0.5*(TMAX - TMIN)*tanh((x)/(L0/10)))
//#define u_BC (Uin*(sq(0.5*L0) - sq(y)))
#define u_BC (0)
u.n[left]  = dirichlet(u_BC);
p[left]    = dirichlet(0);
pf[left]   = dirichlet(0);
f[left]    = dirichlet(1);
fs[left]    = dirichlet(0);
T[left]    = neumann(0);
//T[left]    = dirichlet(Tin + (Tcyl - Tin)*tanh(100*t));
alpha_doc[left] = neumann(0);

u.n[right] = dirichlet(0);
p[right]   = dirichlet(0);
pf[right]  = dirichlet(0);
T[right]    = neumann(0);
f[right]   = dirichlet(1);
fs[right]   = dirichlet(0);
alpha_doc[right] = neumann(0);

u.n[bottom] = dirichlet(0);
p[bottom]   = dirichlet(0);
pf[bottom]  = dirichlet(0);
T[bottom]    = neumann(0);
f[bottom]   = dirichlet(1);
fs[bottom]   = dirichlet(0);
alpha_doc[bottom] = neumann(0);

u.n[top] = dirichlet(0);
p[top]   = dirichlet(0);
pf[top]  = dirichlet(0);
T[top]    = neumann(0);
f[top]   = dirichlet(1);
fs[top]   = dirichlet(0);
alpha_doc[top] = neumann(0);

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
    return 1; // no front, no bubbles
}

double xmin_center = 0, xmax_center = 0;


double geometry(double x, double y, double z){
    coord mypoint = {x,y,z};
    coord pnt_dist;
    double phi = -HUGE;
    foreach_dimension() pnt_dist.x = mypoint.x;
    phi = max(phi, (1 - mynorm(pnt_dist)));
    return phi;
}

void solid_func(scalar fs, face vector fs_face){
    vertex scalar phi[];
    foreach_vertex() {
        phi[] = geometry(x, y, z);
    }
    boundary ({phi});
    fractions (phi, fs, fs_face);
}

event init (t = 0) {
    char name[300];
    sprintf (name, "restart_%s", subname);

    if (!restore (file = name)) {
        fprintf(ferr, "The file %s can not be successfully read! Iniitial conditions are set\n", name);
        int it = 0;
        scalar f_smoothed[], fs_smoothed[];
        do {
            solid_func(fs, fs_face);
            fraction(f, bubbles(x, y, z));
            filter_scalar(f, f_smoothed);
            filter_scalar(fs, fs_smoothed);
        }while (adapt_wavelet({f_smoothed, fs_smoothed}, (double []){feps, feps}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
        foreach() {
            T[] = var_hom(f[], fs[], Tin, Tam, Tcyl);
            alpha_doc[] = 0;
            u.x[] = u_BC*(1 - fs[]); //u_BC*f[];// 0; // penalization will work
            u.y[] = 0;
            un[] = u.x[];
            if (rho1 || rho2){
                rhov[] = var_hom(f[], fs[], rho1, rho2, rho3);
            }
        }
        boundary ({T, alpha_doc, u, rho});
        foreach_face(){
            double ff1 = (f[] + f[-1])/2.;
            double ff2 = (fs[] + fs[-1])/2.; //solid
            if (kappa1 || kappa2)
                kappav.x[] = var_hom(ff1, ff2, kappa1, kappa2, kappa3);
            if (mu1 || mu2) {
                double Tf = (T[] + T[-1])/2.;
                double alphaf = (alpha_doc[] + alpha_doc[-1])/2.;
                face vector muv = mu;
                muv.x[] = fm.x[] * mu(ff1, ff2, alphaf, Tf); //((1 - clamp(fs,0.,1.)) < SEPS);//
            }
        }
		boundary((scalar *){kappa, mu});
        event("vtk_file");
    }else{
        FILE *popen(const char *cmd_str, const char *mode);
        int pclose(FILE *stream);
        FILE *cmd;
        char result[5000];
        char cmd_str[200];
        strcpy(cmd_str, "grep \"vtk: iter_fp\" ");
        strcat(cmd_str, logname);
        strcat(cmd_str, " | awk \'{print $7}\' ");
        fprintf(ferr, "grep timesteps: %s", cmd_str);
        cmd = popen(cmd_str, "r");
        if (cmd == NULL) {
            fprintf(ferr, "Error in opening log file and searching timesteps");
            perror("popen");
            exit(EXIT_FAILURE);
        }
        int k = 0;
        while (fgets(result, sizeof(result), cmd)) {
            printf ("%s", result);
            file_timesteps[k++] = atof(result);
        }
        cmd_str[0] = 0;
        strcpy(cmd_str, "grep \"vtk: iter_fp\" ");
        strcat(cmd_str, logname);
        strcat(cmd_str, " | tail -1 | awk \'{print $4}\' ");
        fprintf(ferr, "grep iter_fp: %s", cmd_str);
        cmd = popen(cmd_str, "r");
//        cmd = popen("grep \"vtk: iter_fp\" log | tail -1 | awk \'{print $4}\'", "r");
        if (cmd == NULL) {
            fprintf(ferr, "Error in opening log file and searching iter_fp");
            perror("popen");
            exit(EXIT_FAILURE);
        }
        fgets(result, sizeof(result), cmd);
        iter_fp = atoi(result) + 1;
        fprintf(ferr, "Read iter_fp+1: %d\n", iter_fp);
        pclose(cmd);
    }
}

event properties(i+=100){

    double mu1_tmp, alpha_doc_max=-1e+10, T_min=1e+10;

    foreach(reduction(max:alpha_doc_max) reduction(min:T_min)){
        if (T[] < T_min) T_min = T[];
        if (alpha_doc[] > alpha_doc_max) alpha_doc_max = alpha_doc[];
    }

#if REACTION_MODEL != NO_REACTION_MODEL
    mu1_tmp = mu(1, 0, alpha_doc_max, T_min); //inside resin
#else
    mu1_tmp = mu(1, 0, 0, T_min); //inside resin
#endif
    mu3 = mu1_tmp*MuRS;
    fprintf(ferr, "mu3=%g was updated\n", mu3);
}

event set_dtmax (i++) {
	  RELATIVE_RES_TOLERANCE = 0.01; //fabs(res1 - res2)/(res1 + 1e-30) < RELATIVE_RES_TOLERANCE
    DT *= 1.05;
    DT = min(DT, maxDT);
//    fprintf(ferr, "TIMEMAX set_dtmax: tnext= %g, t=%g, DT=%g, dt=%g\n", tnext, t, DT, dt);
}


/**
The gravity vector is aligned with the channel and viscosity is
unity. */

//event acceleration (i++) {
//  if (fabs(Ggrav_ndim) > 1e-10){
//  	foreach_face(x)	av.x[] = Ggrav_ndim;
//  }
//}

event advection_term(i++){
    TOLERANCE = TOLERANCE_P;
}

event viscous_term(i++){
    TOLERANCE = TOLERANCE_V;
    m_bp = 1;//max(3 - 0.1*i, 0.4);
    m_bp_T = m_bp;
    double mindelta = L0/pow(2, maxlevel);
    double nu_bp = mu1/rho1;
    eta_s = sq(m_bp*mindelta)/nu_bp;
    eta_T = sq(m_bp_T*mindelta)/chi_conductivity;
//    fprintf(ferr, "VISC: m=%g, mindelta=%g, nu=%g, chi_conductivity=%g, eta_s=%15.12g, eta_T=%15.12g\n", m_bp, mindelta, nu_bp, chi_conductivity, eta_s, eta_T);
}

event projection(i++){
    TOLERANCE = TOLERANCE_P;
}

event end_timestep(i++){
	TOLERANCE = TOLERANCE_T;
}

FILE *fp;
int firstWrite = 0;
event logoutput(t+=0.01){
        char name_vtu[1000];
        coord loc[1];
        double v[2];
        loc[0].x = rslice;
        loc[0].y = 0;
        loc[0].z = 0;
        sprintf(name_vtu, "cylinder_polymerization_basilisk_r=%g.csv", rslice);

        if (firstWrite == 0){
            fp = fopen(name_vtu, "w");
            if(pid() == 0) fprintf (fp, "t,T,alpha\n");
            fclose(fp);
            firstWrite++;
        }
        fp = fopen(name_vtu, "a");
        interpolate_array ((scalar*) {T, alpha_doc}, loc, 1, v, true);
        if(pid() == 0)  fprintf (fp, "%g,%g,%g\n", t, v[0], v[1]);
        fclose(fp);
}

event vtk_file (t += dt_vtk)
//event vtk_file (i += 1)
{
    char name[300];
    sprintf (name, "vtk_%s", subname);
    scalar l[];
    calc_scalar_from_face(mu, mu_cell);
    foreach() {
        l[] = level;
    }
    boundary((scalar *){l});
    fprintf(ferr, "output_vtu_MPI\n");
#ifdef DEBUG_BRINKMAN_PENALIZATION
    output_vtu_MPI(name, (iter_fp) ? t + dt : 0, list = (scalar *) {T, alpha_doc, p, fs, f, l, rhov, mu_cell, my_kappa},
    vlist = (vector *) {u, g, uf, a, dbp, total_rhs, residual_of_u, divtauu, fs_face, alpha});
#else
    output_vtu_MPI(name, (iter_fp) ? t + dt : 0, list = (scalar *) {T, alpha_doc, fs, l, rhov, mu_cell});
#endif
}



/**
We adapt according to the error on the embedded geometry, velocity and
tracer fields. */

#define ADAPT_SCALARS {rhov, fs, T, alpha_doc, mu_cell}
#define ADAPT_EPS_SCALARS {rhoeps, fseps, Teps, aeps, mueps}
event adapt (i++){
  calc_scalar_from_face(mu, mu_cell);
	double eps_arr[] = ADAPT_EPS_SCALARS;
	MinMaxValues((scalar *) ADAPT_SCALARS, eps_arr);
	adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
    count_cells(t, i);
}

//event snapshot (i += snapshot_i)
//{
//    char name[300];
//    if (t==0)
//        sprintf (name, "dump_%s-0.0", subname);
//    else
//        sprintf (name, "dump_%s-%g", subname, t);
//    dump (file = name);
//}

event check_fail(i += 20){
    double umax = 0;
    foreach(){
        if ((u.x[] != u.x[]) || (umax > 10e+10)) {
            fprintf(ferr, "NAN values, u.x=%g, umax=%g", u.x[], umax);
            return 9;
        }
    }
}
//event stop(i=100);
event stop(t = 10*L0 / Uin);


//FILE *fp;
//event logoutput(t+=0.01){
//    int Np = 100;
//    fprintf(ferr, "logoutput:");
//    coord * loc = malloc (Np*sizeof(coord));
//    double * v = malloc (2*Np*sizeof(double));
//
//    char name_vtu[1000];
//    sprintf(name_vtu, "alpha_T_%s.csv", subname);
//    if (i == 0){
//        fp = fopen(name_vtu, "w");
//        if(pid() == 0)
//            fprintf (fp, "t,x,alpha_doc,T\n");
//    }
//
//    for (int ii = 0; ii < Np; ii++){
//        loc[ii].x = X0 + L0*ii/(Np-1);
//        loc[ii].y = 0;
//    }
//    interpolate_array ((scalar*) {alpha_doc, T}, loc, Np, v, true);
//    for (int ii = 0; ii < Np; ii++){
//        if(pid() == 0)
//            fprintf (fp, "%g,%g,%g,%g\n", t, loc[ii].x, v[2*ii], v[2*ii+1]);
//    }
//}