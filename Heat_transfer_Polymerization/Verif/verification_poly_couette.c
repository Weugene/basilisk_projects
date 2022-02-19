#define BRINKMAN_PENALIZATION 1
#define DEBUG_MINMAXVALUES
//#define DEBUG_BRINKMAN_PENALIZATION
//#define DEBUG_MODE_POISSON
//#define DEBUG_HEAT
#define REACTION_MODEL REACTION_MODEL_NON_AUTOCATALYTIC
//#define DEBUG_OUTPUT_VTU_MPI
#define FILTERED
#define JACOBI 1
#define CORRECT_UF_FLUXES 1
#define STICKY_SOLID 1
#define RELATIVE_RES
#define EPS_MAXA 2 // adapt based on Max-Min value
//#define PRINT_ALL_VALUES
//#define STOKES
#define T_DIRICHLET_BC 1

#define mu(f, fs, alpha_doc, T) (mu2 + (muf1(alpha_doc, T) - mu2)*clamp(f,0.,1.))

scalar mu_cell[];
scalar Phi_visc[], Phi_src[];
face vector fs_face[];
vector target_Uv[];
scalar T_targetv[];
#include "centered-weugene.h"
#include "rheology_model.h"
#include "utils-weugene.h"
#include "output_vtu_foreach.h"

int snapshot_i = 1000;
double dt_vtk = 0.01;
double Uin, Tin, T_solid, Tam;
double RhoR, RhoRS, MuR, MuRS, CpR, CpRS, KappaR, KappaRS;
double Rho1, Rho2, Rho3;
double Mu0, Mu1, Mu2, Mu3;
double Kappa1, Kappa2, Kappa3;
double CP1, CP2, CP3;
double Ggrav, Ggrav_ndim;
double channel_diam, channel_length;
double Re; //Reynolds
double Pe; //Peclet number
double Pr; //Prandtl number
double Fr; //Froude number Fr = sqrt(u^2/(g*cyl_diam))
double G;
double Umean;
double layer_velocity, layer_heat;
double x_init = 2, Dx_min, dx_min;
int maxlevel = 9;
int minlevel = 5;
int LEVEL = 7;
double maxDT, maxDT0;
double mindelta, mindelta0;
double mu_max = 0, nu_max = 0;
int adapt_method = 1; // 0 - traditional, 1 - using limitation, 2 - using array for maxlevel
double m_hump = 3;
double feps = 1e-10, fseps = 1e-10, ueps = 1e-6, rhoeps = 1e-10, Teps = 1e-6, aeps = 1e-6, mueps=1e-6;
double TOLERANCE_P = 1e-7, TOLERANCE_V = 1e-8, TOLERANCE_T = 1e-7;

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
    DT = 1e-6;
    maxDT0 = 2.5e-3;

    N_smooth = 1; //three-phase-rheology.h
#ifdef STOKES
    stokes = true;
    stokes_heat = false;
#endif
    m_bp = 2; // Brinkman layer resolution
    m_bp_T = 2;
    fprintf(ferr, "Stokes flow:stokes=%d, stokes_heat=%d, uf.x correction. Solid kappa, Cp are the same as fluid 1\n", stokes, stokes_heat);
#if T_DIRICHLET_BC != 0
    fprintf(ferr, "T_DIRICHLET_BC: 1\n");
#endif
    sprintf(subname, "saturated_1");
    sprintf(logname, "log_saturated_1");
// Physical parameters
	Uin = 1e-2;
    channel_diam = 15e-6;
	Tin = 300, T_solid = 400, Tam = 300;
//	Gorthala, R., Roux, J. A., & Vaughan, J. G. (1994). Resin Flow, Cure and Heat Transfer Analysis for Pultrusion Process. Journal of Composite Materials, 28(6), 486–506. doi:10.1177/002199839402800601  sci-hub.se/10.1177/002199839402800601
//EPON 9310
	Htr = 355000; //J/kg
	Arrhenius_const = 80600;//1/s //originally in Gothala: 80600
	Ea_by_R = 12000/8.314;// Kelvin //origianally  64000/8.314;
    n_degree = 1.2;

	Eeta_by_Rg = 3.76e+4/8.314;// Kelvin Epon 9310,Safonov page 21
	chi = 20;
    Rho1 = 1200, Rho2 = 1.092, Rho3 = 1200;//1790;//air at 23 C Graphite
    CP1 = 1255, CP2 = 1006, CP3 = 1255;//712;//J/(kg*K) !!!changed
    Kappa1 = 0.002, Kappa2 = 0.02535, Kappa3 = 0.002;// 0.2 --- 8.70;//W/(m*K) !!!changed
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

    channel_length = (2.0/(1.0 - pow(2.0, 1 - maxlevel)))*channel_diam;
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
    Pe = CP1*Rho1*Uin*channel_diam/Kappa1;
    Pr = Mu1*CP1/Kappa1;
// Dimensionless parameters are chosen channel_diam, rho1, Cp1, Tin, Uin
    size(channel_length/channel_diam);
    mindelta = L0/pow(2, maxlevel);
    mindelta0 = L0/pow(2, 8);
    maxDT = maxDT0*sq(mindelta/mindelta0);  // h^2/DT = h0^2/DT0
	RhoR = Rho2/Rho1, RhoRS = Rho3/Rho1;
	MuR = Mu2/Mu1, MuRS = 1;// Mu3/Mu1; //
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
    layer_velocity = channel_diam/sqrt(Re);
    layer_heat = channel_diam/sqrt(Pe);
    eta_s = sq(m_bp*mindelta)/(mu1/rho1);
    eta_T = sq(m_bp_T*mindelta)/chi_conductivity;
    init_grid(1 << maxlevel);
    origin (-L0/2, -L0/2.);
    periodic(right);
//    periodic(top);
	fprintf(ferr,"Dim-less vars: mu0=%g, mu1=%g, mu2=%g, mu3=%g,\n"
                 "               rho1=%g, rho2=%g, rho3=%g,\n"
                 "               nu1=%g, nu2=%g, nu3=%g,\n"
				 "               kappa1=%g, kappa2=%g, kappa3=%g,\n"
                 "               Cp1=%g, Cp2=%g, Cp3=%g,\n"
                 "               RhoR=%g, RhoRS=%g, MuR=%g, MuRS=%g,\n"
                 "               KappaR=%g, KappaRS=%g, CpR=%g, CpRS=%g\n"
				 "               Uin=%g, Tin=%g, T_solid=%g, Tam=%g,\n"
				 "               Htr=%g, Arrhenius_const=%g, Ea_by_R=%g, Eeta_by_Rg=%g,\n"
				 "               L0=%g, channel_diam=%g, maxDT0=%g for maxlevel=8, maxDT=%g for maxlevel=%d\n"
				 "               Ggrav_ndim=%g Uin=%g,\n"
                 "               layer_velocity=%g layer_heat=%g\n"
                 "Dim-less nums: Re=%g,  Fr=%g, Pe=%g, Pr=%g\n"
                 "Solver:        DTmax=%g, CFL=%g, CFL_ARR=%g, NITERMIN=%d,  NITERMAX=%d,\n"
                 "               TOLERANCE_P=%g, TOLERANCE_V=%g, TOLERANCE_T=%g\n"
                 "ADAPT:         minlevel=%d,  maxlevel=%d, feps=%g, fseps=%g, ueps=%g, Teps=%g, aeps=%g\n"
                 "OUTPUT:        dt_vtk=%g,    number of procs=%d\n",
				mu0, mu1, mu2, mu3, rho1, rho2, rho3,
				mu1/rho1, mu2/rho2, mu3/rho3,
				kappa1, kappa2, kappa3,
				Cp1, Cp2, Cp3,
                RhoR, RhoRS, MuR, MuRS,
                KappaR, KappaRS, CpR, CpRS,
				Uin, Tin, T_solid, Tam,
				Htr, Arrhenius_const, Ea_by_R, Eeta_by_Rg,
				L0, channel_diam, maxDT0, maxDT, maxlevel,
				Ggrav_ndim, Uin,
                layer_velocity, layer_heat,
                Re, Fr, Pe, Pr,
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

#define hump(x, xh) ((fabs(x) < xh) ? (1.0/pow(xh, 6))*cube(sq(xh) - sq(x)) : 0.0)
// #define u_BC (0) // steady-state
//#define u_BC (min(max(Uin*(y + 0.5), 0), Uin))  // limited linear function
#define u_BC ( (fabs(y) <= 0.5) ? 0.5 * Uin - 0.5 * Uin * cos(pi*(y + 0.5)) :  \
               (y < -0.5) ? 0 : Uin   ) // Smoothed
//#define T_BC (min(max((T_solid - Tin)*y + 0.5*(T_solid + Tin), Tin), T_solid)) // limited linear function
//#define T_BC min(T_solid + (T_solid - Tin)*tanh( 100*(y - 0.5) ), T_solid)
//#define T_BC ( (fabs(y) <= 0.5) ? 0.5 * (T_solid + Tin) - 0.5 * (T_solid - Tin) * cos(2*pi*y) :  T_solid   )
//#define T_BC ( (fabs(y) <= 0.5) ? 0.5 * (T_solid + Tin) - 0.5 * (T_solid - Tin) * cos(pi*(y + 0.5)) :  \
//               (y < -0.5) ? Tin : T_solid   ) // Smoothed
//#define T_BC ( (fabs(y) <= 0.5) ? 0.5 * (T_solid + Tin) + 0.5 * (T_solid - Tin) * cos(2*pi*m_hump*(y + 0.5)) :  \
//               (y < -0.5) ? T_solid : T_solid   ) // Smoothed with hump
#define T_BC ( T_solid + (Tin - T_solid)*(hump(y - 0.05, 0.05) + hump(y + 0.05, 0.05) )) // Smoothed with 2 humps


u.n[bottom] = dirichlet(0);
u.t[bottom] = dirichlet(u_BC);
p[bottom]   = neumann(0);
pf[bottom]  = neumann(0);
T[bottom]    = dirichlet(T_BC);
f[bottom]   = dirichlet(1);
fs[bottom]   = dirichlet(1);
alpha_doc[bottom] = neumann(0);

u.n[top] = dirichlet(0);
u.t[top] = dirichlet(u_BC);
p[top]   = neumann(0);
pf[top]  = neumann(0);
T[top]    = dirichlet(T_BC);
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
        target_Uv.x[] = ( y > 0 ) ? Uin : 0;
        target_Uv.y[] = 0;
        T_targetv[] = T_solid;
    }
    boundary((scalar *){target_Uv, T_targetv});
}

double exact_solution(double x, double y, double z, double t){
    double se = sqrt(eta_T);
    double D1 = (T_solid - Tin)/(1.0 + 2.0*se);
    double D2 = 0.5*(Tin + T_solid);
    double A2 = se*D1;
    double B1 = -se*D1;

    if (y < -0.5)
        return A2*exp((y + 0.5)/se) + Tin;
    else if (y < 0.5)
        return D1*y + D2;
    else
        return B1*exp(-(y - 0.5)/se) + T_solid;
}

event init (t = 0) {
    char name[300];
    sprintf (name, "restart_%s", subname);

    if (!restore (file = name)) {
        fprintf(ferr, "The file %s can not be successfully read! Initial conditions are set\n", name);
         int it = 0;
         scalar fs_smoothed[];
         do {
             solid_func(fs, fs_face, target_Uv, T_targetv);
             filter_scalar(fs, fs_smoothed);
             foreach() {
                 T[] = T_BC;
                 u.x[] = u_BC;
             }
             boundary ((scalar *){T, u.x});
         }while (adapt_wavelet({fs_smoothed, T, u.x}, (double []){feps, Teps, ueps}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
        solid_func(fs, fs_face, target_Uv, T_targetv);
        foreach() {
            f[] = 1;
            alpha_doc[] = 0;
            u.y[] = 0;
            if (rho1 || rho2){
                rhov[] = var_hom(f[], fs[], rho1, rho2, rho3);
            }
        }
        boundary ((scalar *){f, alpha_doc, u.y, rhov});
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
        muv[] = mu(f[], fs[], alpha_doc[], T[]);
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

event set_dtmax (i++) {
    RELATIVE_RES_TOLERANCE = 0.01; //fabs(res1 - res2)/(res1 + 1e-30) < RELATIVE_RES_TOLERANCE
    DT *= 1.05;
    DT = min(DT, maxDT);
    fprintf(ferr, "TIMEMAX i=%d set_dtmax: tnext= %g, t=%g, DT=%g, dt=%g\n", i, tnext, t, DT, dt);
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
double yslice;
const int Nvar = 4;
event logoutput(t += 0.01){
    yslice = 0.4;
//    int Nh = (yslice - 0.5*mindelta)/mindelta;
//    yslice = Nh*mindelta + 0.5*mindelta;
    char name_vtu[1000];
    coord loc[1];
    double v[Nvar];
    loc[0].x = 0;
    loc[0].y = yslice;
    loc[0].z = 0;
    sprintf(name_vtu, "couette_polymerization_basilisk_x=%g.csv", yslice);

    if (firstWrite == 0){
        fp = fopen(name_vtu, "w");
        if(pid() == 0) fprintf (fp, "t,T,alpha,u,mu,Te\n");
        fclose(fp);
        firstWrite++;
    }
    fp = fopen(name_vtu, "a");
    interpolate_array ((scalar*) {T, alpha_doc, u.x, mu_cell}, loc, 1, v, true);
    if(pid() == 0)  fprintf (fp, "%g,%g,%g,%g,%g,%g\n", t, v[0], v[1], v[2],v[3], exact_solution(0, yslice, 0, t));
    fclose(fp);
}

const int Ninterp = 1001;
event logoutput2(t += 0.01){
    if (t < 2){
        char name_vtu[1000];
        coord loc[Ninterp];
        double v[Nvar*Ninterp];
        for (int i = 0; i < Ninterp; i++){
            loc[i].y = -0.5 + i/(Ninterp - 1.0);
            loc[i].x = loc[i].z = 0;
        }
        sprintf(name_vtu, "couette_polymerization_basilisk_tfix=%g.csv", t);
        fp = fopen(name_vtu, "w");
        if (pid() == 0) fprintf (fp, "x,T,alpha,u,mu,Te\n");
        interpolate_array ((scalar*) {T, alpha_doc, u.x,mu_cell}, loc, Ninterp, v, true);
        if (pid() == 0){
            for (int i = 0; i < Ninterp; i++){
                int ii = i*Nvar;
                fprintf (fp, "%g,%g,%g,%g,%g,%g\n", loc[i].y, v[ii], v[ii+1], v[ii+2], v[ii+3], exact_solution(0, loc[i].y, 0, t));
            }
        }
        fclose(fp);
    }
}

void calcPhiVisc (vector u, face vector uf, scalar T, scalar alpha_doc, scalar f, scalar fs, scalar Phi_visc, scalar Phi_src){
    foreach(){
        double grad2 = 0.0;
        foreach_dimension()
        grad2 += sq((uf.x[1] - uf.x[]));
        Phi_visc[] = 2*grad2/sq(Delta)
#if dimension >= 2
            + sq(0.5*(u.x[0,1] - u.x[0,-1])/Delta) // dv_x/dy
                + sq(0.5*(u.y[1] - u.y[-1])/Delta)  // dv_y/dx
#endif
#if dimension > 2
            + sq(0.5*(u.z[0,1] - u.z[0,-1])/Delta)// dv_z/dy
                + sq(0.5*(u.y[0,0,1] - u.y[0,0,-1])/Delta) // dv_y/dz
                + sq(0.5*(u.x[0,0,1] - u.x[0,0,-1])/Delta) // dv_x/dz
                + sq(0.5*(u.z[1] - u.z[-1])/Delta) // dv_z/dx
#endif
                ;
        Phi_visc[] *= mu(f[], fs[], alpha_doc[], T[]);
        Phi_src[] = f[] * KT(T[]) * FR(alpha_doc[]);
    }
    boundary((scalar *){Phi_visc, Phi_src});
}

event vtk_file (t += dt_vtk)
//event vtk_file (i += 1)
{
    char name[300];
    sprintf (name, "vtk_%s", subname);
    scalar l[], Te[];
    foreach() {
        l[] = level;
        Te[] = exact_solution(x,y,z,t);
    }
    boundary((scalar *){l, Te});
    calcPhiVisc (u, uf, T, alpha_doc, f, fs, Phi_visc, Phi_src);
#ifdef DEBUG_BRINKMAN_PENALIZATION
    output_vtu_MPI(name, (iter_fp) ? t + dt : 0, list = (scalar *) {T, T_targetv, target_Uv, alpha_doc, p, fs, f, l, rhov, mu_cell},
    vlist = (vector *) {u, g, uf, dbp, total_rhs, residual_of_u, divtauu, fs_face, alpha});
#else
    output_vtu_MPI(name, (iter_fp) ? t + dt : 0, list = (scalar *) {T, alpha_doc, fs, l, mu_cell, Phi_visc, Phi_src, Te}, vlist = (vector *){u});
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
    if (i%100) count_cells(t, i);
}

event snapshot (t += 0.5)
{
    char name[300];
    if (t==0)
        sprintf (name, "dump_%s-0.0", subname);
    else
        sprintf (name, "dump_%s-%g", subname, t);
    dump (file = name);
}

event check_fail(i += 100){
    double umax = 0;
    foreach(){
        if ((u.x[] != u.x[]) || (umax > 10e+10)) {
            fprintf(ferr, "NAN values, u.x=%g, umax=%g", u.x[], umax);
            return 9;
        }
    }
}

event stop(t = 2);
