#define BRINKMAN_PENALIZATION 1
//#define DEBUG_MINMAXVALUES
//#define DEBUG_BRINKMAN_PENALIZATION
#define DEBUG_MODE_POISSON
//#define DEBUG_HEAT
#define REACTION_MODEL REACTION_MODEL_NON_AUTOCATALYTIC
//#define DEBUG_OUTPUT_VTU_MPI
#define FILTERED
#define JACOBI 1
//#define CORRECT_UF_FLUXES 1
//#define STICKY_SOLID 1
#define RELATIVE_RES
#define EPS_MAXA 2 // adapt based on Max-Min value
//#define PRINT_ALL_VALUES
#define STOKES
#define T_DIRICHLET_BC 1

scalar mu_cell[];
face vector fs_face[];
vector target_Uv[];
scalar T_targetv[];
#include "centered-weugene.h"
#include "rheology_model.h"
#include "utils-weugene.h"
#include "output_vtu_foreach.h"

int snapshot_i = 1000;
double dt_vtk = 0.1;
double Uin, Tin, T_solid, Tam;
double RhoR, RhoRS, MuR, MuRS, CpR, CpRS, KappaR, KappaRS;
double Rho1, Rho2, Rho3;
double Mu0, Mu1, Mu2, Mu3;
double Kappa1, Kappa2, Kappa3;
double CP1, CP2, CP3;
double Ggrav, Ggrav_ndim;
double channel_diam, channel_length;
double Re; //Reynolds
double Fr; //Froude number Fr = sqrt(u^2/(g*v))
double G;
double Umean, PLeft=10;
double x_init = 2, Dx_min, dx_min;
int maxlevel = 9;
int minlevel = 5;
int LEVEL = 7;
double maxDT, maxDT0;
double mu_max = 0, nu_max = 0;
int adapt_method = 1; // 0 - traditional, 1 - using limitation, 2 - using array for maxlevel
double feps = 1e-10, fseps = 1e-10, ueps = 1e-3, rhoeps = 1e-10, Teps = 1e-4, aeps = 1e-3, mueps=1e-3;
double TOLERANCE_P = 1e-7, TOLERANCE_V = 1e-8, TOLERANCE_T = 1e-7;
double mindelta, mindelta0;

char subname[150], logname[200];
int main(int argc, char * argv[]) {
    fprintf(ferr,"correction of uf\n");
    fprintf(ferr, "./a.out T_solid, Tin, maxlevel, iter_fp, "
                  "TOLERANCE_P, TOLERANCE_V, TOLERANCE_T, Htr, Arrhenius_const, Ea_by_R, subname\n");
    TOLERANCE = 1e-7;
    NITERMIN = 1;
    NITERMAX = 100;
    CFL = 0.4;
    CFL_ARR = 0.5;
    DT = 1e-5;
    maxDT0 = 1e-2;

    N_smooth = 1; //three-phase-rheology.h
#ifdef STOKES
    stokes = true;
    stokes_heat = true;
#endif
    m_bp = 0.5; // Brinkman layer resolution
    m_bp_T = 0.5;
    fprintf(ferr, "Stokes flow:stokes=%d, stokes_heat=%d, uf.x correction. Rho2 is 10x and Mu2 is 100x higher\n", stokes, stokes_heat);
#if T_DIRICHLET_BC != 0
    fprintf(ferr, "T_DIRICHLET_BC: 1\n");
#endif
    sprintf(subname, "saturated_1");
    sprintf(logname, "log_saturated_1");
// Physical parameters
	Uin = 1e-2;
    channel_diam = 15e-6, channel_length = 10*channel_diam;
	Tin = 300, T_solid = 400.15, Tam = 300;
//	Gorthala, R., Roux, J. A., & Vaughan, J. G. (1994). Resin Flow, Cure and Heat Transfer Analysis for Pultrusion Process. Journal of Composite Materials, 28(6), 486–506. doi:10.1177/002199839402800601  sci-hub.se/10.1177/002199839402800601
//EPON 9310
	Htr = 0;//355000; //J/kg
	Arrhenius_const = 80600;//1/s //originally in Gothala: 80600
	Ea_by_R = 12000/8.314;// Kelvin //origianally  64000/8.314;
    n_degree = 1.2;

	Eeta_by_Rg = 3.76e+4/8.314;// Kelvin Epon 9310,Safonov page 21
	chi = 20;
    Rho1 = 1200, Rho2 = 1.092, Rho3 = 1790;//air at 23 C Graphite
    CP1 = 1255, CP2 = 1006, CP3 = 712;//J/(kg*K)
    Kappa1 = 0.2, Kappa2 = 0.02535, Kappa3 = 8.70;//W/(m*K)
    Ggrav = 0; // m/s^2

    if (argc > 1)
        T_solid = atof(argv[1]);
    if (argc > 2)
        Tin = atof(argv[2]);
    Mu0 = 3.85e-7, Mu1 = Mu0*exp(Eeta_by_Rg/Tin), Mu2 = 1.963e-5, Mu3 = Mu0*exp(Eeta_by_Rg/T_solid)*Rho3/Rho1;//air at 50C //Mu2 = 1.963e-5
	if (argc > 3)
        maxlevel = atoi(argv[3]);
    if (argc > 4)
        iter_fp = atoi(argv[4]);
	if (argc > 5)
        TOLERANCE_P = atof(argv[5]);
    if (argc > 6)
        TOLERANCE_V = atof(argv[6]);
    if (argc > 7)
        TOLERANCE_T = atof(argv[7]);
	if (argc > 8)
		Htr = atof(argv[8]);
	if (argc > 9)
		Arrhenius_const = atof(argv[9]);
	if (argc > 10)
		Ea_by_R = atof(argv[10]);
    if (argc > 11) {
        strcpy(subname, argv[11]);
        sprintf (logname, "log%s", subname);
    }

	fprintf(ferr,
                 "Props(SI): Mu0=%g, Mu1=%g, Mu2=%g, Mu3=%g, Rho1=%g, Rho2=%g,  Rho3=%g,\n"
                 "           nu1=%g, nu2=%g, nu3=%g,\n"
                 "           Kappa1=%g, Kappa2=%g, Kappa3=%g, CP1=%g, CP2=%g, CP3=%g,\n"
				 "           Uin=%g, time*=%g, Tin=%g, T_solid=%g, Tam=%g\n"
                 "           Htr=%g, Arrenius=%g, Ea_by_R=%g, n_deg=%g, m_deg=%g\n"
				 "           Eeta_by_Rg=%g, chi=%g\n"
                 "Geometry: channel_diam=%g,  domainSize=%g\n",
                 Mu0, Mu1, Mu2, Mu3, Rho1, Rho2, Rho3,
                 Mu1/Rho1, Mu2/Rho2, Mu3/Rho3,
                 Kappa1, Kappa2, Kappa3, CP1, CP2, CP3,
                 Uin, channel_diam/Uin, Tin, T_solid, Tam,
                 Htr, Arrhenius_const, Ea_by_R, n_degree, m_degree,
                 Eeta_by_Rg, chi,
                 channel_diam, channel_length);
// Dimensionless numbers
	Re = Uin*channel_diam*Rho1/Mu1;
	Fr = sqrt(sq(Uin)/(Ggrav*channel_diam + SEPS));
// Dimensionless parameters are chosen channel_diam, rho1, Cp1, Tin, Uin
    size(channel_length/channel_diam);
    mindelta = L0/pow(2, maxlevel);
    mindelta0 = L0/pow(2, 8);
    maxDT = maxDT0*sq(mindelta/mindelta0);  // h^2/DT = h0^2/DT0
	RhoR = Rho2/Rho1, RhoRS = Rho3/Rho1;
	MuR = Mu2/Mu1, MuRS = Mu3/Mu1;// Mu3/Mu1; // Mu3/Mu1 = Rho3/Rho1 => RhoRS leads bad result?
	CpR = CP2/CP1, CpRS = CP3/CP1;
	KappaR = Kappa2/Kappa1, KappaRS = Kappa3/Kappa1;
    rho1 = 1; rho2 = RhoR; rho3 = RhoRS;
    mu0 = (1./Re)*(Mu0/Mu1); mu1 = (1./Re); mu2 = mu1*MuR; mu3 = mu1*MuRS; mu_eff = 0;
//    mu3 = mu1*(rho3/rho1)*sq(1.0/(m_bp*dx_min));
	Cp1 = 1; Cp2 = Cp1*CpR; Cp3 = Cp1*CpRS;//J/(kg*K)
	kappa1 = Kappa1/(Rho1*CP1*channel_diam*Uin + SEPS), kappa2 = kappa1*KappaR, kappa3 = kappa1*KappaRS;
    chi_conductivity = kappa1 / (rho1 * Cp1);
    Htr /= CP1*Tin;
	Arrhenius_const *= channel_diam/(Uin + SEPS);
	Ea_by_R /= Tin;
	Eeta_by_Rg /= Tin;
	Ggrav_ndim = 1./sq(Fr);
	Uin = 1;
    Tam /= Tin;
	T_solid /= Tin;
	Tin  /= Tin;
	channel_diam = 1;
    init_grid(1 << maxlevel);
    origin (-L0/2, -L0/2.);
    periodic(right);
//    periodic(top);
	fprintf(ferr,"Dim-less vars: mu0=%g, mu1=%g, mu2=%g, mu3=%g,\n"
                 "               rho1=%g, rho2=%g, rho3=%g,\n"
                 "               nu1=%g, nu2=%g, nu3=%g,\n"
				 "               kappa1=%g, kappa2=%g, kappa3=%g,\n"
                 "               Cp1=%g, Cp2=%g, Cp3=%g,\n"
				 "               Uin=%g, Tin=%g, T_solid=%g, Tam=%g,\n"
				 "               Htr=%g, Arrhenius_const=%g, Ea_by_R=%g, Eeta_by_Rg=%g,\n"
				 "               L0=%g, channel_diam=%g, maxDT0=%g for maxlevel=8, maxDT=%g for maxlevel=%d\n"
				 "               Ggrav_ndim=%g Uin=%g\n"
                 "Dim-less nums: Re=%g,  Fr=%g\n"
                 "Solver:        DTmax=%g, CFL=%g, CFL_ARR=%g, NITERMIN=%d,  NITERMAX=%d,\n"
                 "               TOLERANCE_P=%g, TOLERANCE_V=%g, TOLERANCE_T=%g\n"
                 "ADAPT:         minlevel=%d,  maxlevel=%d, feps=%g, fseps=%g, ueps=%g, Teps=%g, aeps=%g\n"
                 "OUTPUT:        dt_vtk=%g,    number of procs=%d\n",
				mu0, mu1, mu2, mu3, rho1, rho2, rho3,
				mu1/rho1, mu2/rho2, mu3/rho3,
				kappa1, kappa2, kappa3, Cp1, Cp2, Cp3,
				Uin, Tin, T_solid, Tam,
				Htr, Arrhenius_const, Ea_by_R, Eeta_by_Rg,
				L0, channel_diam, maxDT0, maxDT, maxlevel,
				Ggrav_ndim, Uin,
                Re, Fr,
                DT, CFL, CFL_ARR, NITERMIN, NITERMAX,
                TOLERANCE_P, TOLERANCE_V, TOLERANCE_T,
                minlevel, maxlevel, feps, fseps, ueps, Teps, aeps,
                dt_vtk, npe());

    fs.refine = fs.prolongation = fraction_refine;
    f.refine = f.prolongation = fraction_refine;
    boundary((scalar *){f, fs});
    target_U = target_Uv;
    T_target = T_targetv;
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

    run();
}
// #define u_BC (0) // steady-state
#define u_BC (min(max(Uin*(y + 0.5), 0), Uin))  // limited linear function
//#define T_BC (min(max((T_solid - Tin)*y + 0.5*(T_solid + Tin), Tin), T_solid)) // limited linear function
#define T_BC min(T_solid + (T_solid - Tin)*tanh( 100*(y - 0.5) ), T_solid)

u.n[bottom] = dirichlet(0);
u.t[bottom] = dirichlet(0);
p[bottom]   = neumann(0);
pf[bottom]  = neumann(0);
T[bottom]    = dirichlet(Tin);
f[bottom]   = dirichlet(1);
fs[bottom]   = dirichlet(1);
alpha_doc[bottom] = neumann(0);

u.n[top] = dirichlet(0);
u.t[top] = dirichlet(Uin);
p[top]   = neumann(0);
pf[top]  = neumann(0);
T[top]    = dirichlet(T_solid);
f[top]   = dirichlet(1);
fs[top]   = dirichlet(1);
alpha_doc[top] = neumann(0);

double xmin_center = 0, xmax_center = 0;

double geometry(double x, double y, double z){
    return fabs(y) - 0.5;
}

void solid_func(scalar fs, face vector fs_face, vector target_Uv, scalar T_targetv){
    vertex scalar phi[];
    foreach_vertex() {
        phi[] = geometry(x, y, z);
    }
    boundary ((scalar *){phi});
    fractions (phi, fs, fs_face);

    foreach(){
        target_Uv.x[] = ( y > 0 ) ? u_BC : 0;
        target_Uv.y[] = 0;
        T_targetv[] = ( y > 0 ) ? T_solid : Tin;
    }
    boundary((scalar *){target_Uv, T_targetv});
}


event init (t = 0) {
    char name[300];
    sprintf (name, "restart_%s", subname);

    if (!restore (file = name)) {
        fprintf(ferr, "The file %s can not be successfully read! Iniitial conditions are set\n", name);
         int it = 0;
         scalar fs_smoothed[];
         do {
             solid_func(fs, fs_face, target_Uv, T_targetv);
             filter_scalar(fs, fs_smoothed);
         }while (adapt_wavelet({fs_smoothed}, (double []){feps, feps}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
        solid_func(fs, fs_face, target_Uv, T_targetv);
        foreach() {
            f[] = 1;
            T[] = T_BC;
            T_targetv[] = T_BC;
            alpha_doc[] = 0;
            u.x[] = u_BC; //u_BC*f[];// 0; // penalization will work
            u.y[] = 0;
            if (rho1 || rho2){
                rhov[] = var_hom(f[], fs[], rho1, rho2, rho3);
            }
        }
        boundary ((scalar *){f, T, T_targetv, alpha_doc, u, rhov});
        foreach_face(){
            double ff2 = (fs[] + fs[-1])/2.; //solid
            if (kappa1 || kappa2)
                kappav.x[] = var_hom(1, ff2, kappa1, kappa2, kappa3);
            if (mu1 || mu2) {
                double Tf = (T[] + T[-1])/2.;
                face vector muv = mu;
                muv.x[] = fm.x[] * mu(1, ff2, 0, Tf); //((1 - clamp(fs,0.,1.)) < SEPS);//
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

void calc_mu(scalar muv){
    foreach(){
        double mu1_tmp = muf1(alpha_doc[], T[]);
        muv[] = var_hom(f[], fs[], mu1_tmp, mu2, mu3);
        if( muv[] > mu_max) mu_max = muv[];
    }
    boundary((scalar *){muv});
}

event properties(i++){
    solid_func(fs, fs_face, target_Uv, T_targetv);
    calc_mu(mu_cell);
    mu3 = mu_max*MuRS;
    fprintf(ferr, "mu3=%g was updated, mu_max=%g\n", mu3, mu_max);
}

//event stability(i++){
//    stokes = (i<10) ? true : false;
//}
event set_dtmax (i++) {
    RELATIVE_RES_TOLERANCE = 0.01; //fabs(res1 - res2)/(res1 + 1e-30) < RELATIVE_RES_TOLERANCE
    DT *= 1.05;
    DT = min(DT, maxDT);
     fprintf(ferr, "TIMEMAX set_dtmax: tnext= %g, t=%g, DT=%g, dt=%g\n", tnext, t, DT, dt);
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
    nu_max = mu_max/rho1;
    TOLERANCE = TOLERANCE_V;
    eta_s = sq(m_bp*mindelta)/nu_max;
    eta_T = sq(m_bp_T*mindelta)/chi_conductivity;
   fprintf(ferr, "VISC: m=%g, mindelta=%g, nu=%g, chi_conductivity=%g, eta_s=%15.12g, eta_T=%15.12g\n", m_bp, mindelta, nu_max, chi_conductivity, eta_s, eta_T);
}

event projection(i++){
    TOLERANCE = TOLERANCE_P;
}

event end_timestep(i++){
	TOLERANCE = TOLERANCE_T;
}

FILE *fp;
int firstWrite = 0;
double xslice = 0;
const int Nvar = 4;
event logoutput(t+=0.01){
    char name_vtu[1000];
    coord loc[1];
    double v[Nvar];
    loc[0].x = xslice;
    loc[0].y = 0;
    loc[0].z = 0;
    sprintf(name_vtu, "cylinder_polymerization_basilisk_x=%g.csv", xslice);

    if (firstWrite == 0){
        fp = fopen(name_vtu, "w");
        if(pid() == 0) fprintf (fp, "t,T,alpha,u,mu\n");
        fclose(fp);
        firstWrite++;
    }
    fp = fopen(name_vtu, "a");
    interpolate_array ((scalar*) {T, alpha_doc, u.x, mu_cell}, loc, 1, v, true);
    if(pid() == 0)  fprintf (fp, "%g,%g,%g,%g,%g\n", t, v[0], v[1], v[2],v[3]);
    fclose(fp);
}

double tfix = 0.1;
const int Ninterp = 101;
event logoutput2(t+=1){
    char name_vtu[1000];
    coord loc[Ninterp];
    double v[Nvar*Ninterp];
    for (int i = 0; i < Ninterp; i++){
        loc[i].x = -0.5 + i/(Ninterp - 1.0);
        loc[i].y = loc[i].z = 0;
    }
    sprintf(name_vtu, "cylinder_polymerization_basilisk_tfix=%g.csv", t);
    fp = fopen(name_vtu, "w");
    if (pid() == 0) fprintf (fp, "x,T,alpha,u,mu\n");
    interpolate_array ((scalar*) {T, alpha_doc, u.x,mu_cell}, loc, Ninterp, v, true);
    if (pid() == 0){
        for (int i = 0; i < Ninterp; i++){
            int ii = i*Nvar;
            fprintf (fp, "%g,%g,%g,%g,%g\n", loc[i].x, v[ii], v[ii+1], v[ii+2], v[ii+3]);
        }
    }
    fclose(fp);
}

event vtk_file (t += dt_vtk)
//event vtk_file (i += 1)
{
    char name[300];
    sprintf (name, "vtk_%s", subname);
    scalar l[];
    foreach() {
        l[] = level;
    }
    boundary((scalar *){l});
#ifdef DEBUG_BRINKMAN_PENALIZATION
    output_vtu_MPI(name, (iter_fp) ? t + dt : 0, list = (scalar *) {T, T_targetv, target_Uv, alpha_doc, p, fs, f, l, rhov, mu_cell},
    vlist = (vector *) {u, g, uf, dbp, total_rhs, residual_of_u, divtauu, fs_face, alpha});
#else
    output_vtu_MPI(name, (iter_fp) ? t + dt : 0, list = (scalar *) {T, alpha_doc, fs, l, mu_cell}, vlist = (vector *){u});
#endif
}



/**
We adapt according to the error on the embedded geometry, velocity and
tracer fields. */
#define ADAPT_SCALARS {rhov, T, alpha_doc, u.x}
#define ADAPT_EPS_SCALARS {rhoeps, Teps, aeps, ueps}
event adapt (i++){
	double eps_arr[] = ADAPT_EPS_SCALARS;
	MinMaxValues((scalar *) ADAPT_SCALARS, eps_arr);
	adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
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

event check_fail(i += 100){
    double umax = 0;
    foreach(){
        if ((u.x[] != u.x[]) || (umax > 10e+10)) {
            fprintf(ferr, "NAN values, u.x=%g, umax=%g", u.x[], umax);
            return 9;
        }
    }
}
event stop(t = 3 * L0 / Uin);
//event stop(t = 20);
