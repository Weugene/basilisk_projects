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
int adapt_method = 1; // 0 - traditional, 1 - using limitation, 2 - using array for maxlevel
double feps = 1e-10, fseps = 1e-10, ueps = 1e-2, Teps = 3e-2, aeps = 3e-2;
double TOLERANCE_P = 1e-6, TOLERANCE_V = 1e-7, TOLERANCE_T = 1e-6;
double *R = NULL;
coord *centers = NULL;

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

    fprintf(ferr, "The solver must be run: ./a.out Tcyl, Tin, maxlevel, iter_fp, ratio_Rbmin, ratio_Rbmax, "
                  "ratio_dist_x, ratio_dist_y, ratio_front_x, cyl_x, Nb, Ncx, Ncy, "
                  "TOLERANCE_P, TOLERANCE_V, TOLERANCE_T, Htr, Arrhenius_const, Ea_by_R\n");
    if (argc > 1)
        Tcyl = atof(argv[1]);
    if (argc > 2)
        Tin = atof(argv[2]);
//	Mu1 = 3.85e-7*exp(Eeta_by_Rg/Tin), Mu2 = 1e-4, Mu3 = Rho3*Mu1/Rho1/2.0;
    Mu0 = 3.85e-7, Mu1 = Mu0*exp(Eeta_by_Rg/Tin), Mu2 = 1.963e-5, Mu3 = Mu1*Rho3/Rho1;//air at 50C //Mu2 = 1.963e-5

	if (argc > 3)
        maxlevel = atoi(argv[3]);
    //stokes = true;
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
    if (argc > 9 )
        cyl_x = atof(argv[9]);
	if (argc > 10)
        Nb = atoi(argv[10]);
	if (argc > 11)
        Ncx = atoi(argv[11]);
	if (argc > 12)
        Ncy = atoi(argv[12]);
	if (argc > 13)
        TOLERANCE_P = atof(argv[13]);
    if (argc > 14)
        TOLERANCE_V = atof(argv[14]);
    if (argc > 15)
        TOLERANCE_T = atof(argv[15]);
	if (argc > 16)
		Htr = atof(argv[16]);
	if (argc > 17)
		Arrhenius_const = atof(argv[17]);
	if (argc > 18)
		Ea_by_R = atof(argv[18]);
    if (argc > 19)
        dev_r = atof(argv[19]);
    if (argc > 20)
        develx = atof(argv[20]);
    if (argc > 21)
        devely = atof(argv[21]);
    if (argc > 22)
        shift_x = atof(argv[22]);
    if (argc > 23)
        shift_y = atof(argv[23]);
    if (argc > 24)
        non_saturated = atoi(argv[24]);

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

    int rank, psize, h_len;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    // get rank of this proces
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // get total process number
    MPI_Comm_size(MPI_COMM_WORLD, &psize);
    MPI_Get_processor_name(hostname, &h_len);
    printf("rank:%d size: %d at %s h_len %d\n", rank, psize, hostname, h_len);

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
    return (non_saturated > 0)? front_x - x : 1;// with front

}
double xmin_center = 0, xmax_center = 0;
void calc_centers(coord * centers, double * R){
    int k;
    coord pnt_dist;
    double minvx, maxvx, minvy, maxvy;
    double limMin = cyl_x + cyl_diam;
    int trials = 0;
// generating of Radii and centers of bubbles which are not overlapping
    srand (0); // for consistency
    for (int i = 0; i < Ncx; i++){
        for(int j = 0; j < Ncy; j++){
            k = i*Ncy + j;
            bool flag = 1;
            minvx = limMin + i*dist_x + shift_x*((j%2) > 0), maxvx = minvx + develx;
            minvy = -0.5*L0 + 0.5*dist_y + j*dist_y - shift_y*((i%2) > 0), maxvy = minvy + devely;
//            if ( (maxvx > X0 + L0 - cyl_diam) || (maxvy > Y0 + L0 - 0.5*cyl_diam) || (minvy < Y0 + 0.5*cyl_diam)){
//                // otherwise centers are far away. In such a way I get rid of this points
//                R[k] = 0;
//                foreach_dimension() centers[k].x = 1e+9;
//                continue;
//            }
            trials = 0;
            while(1){
                R[k] = RandMinMax(0.5*cyl_diam, 0.5*cyl_diam + dev_r);
                centers[k].x = RandMinMax(minvx, maxvx);
                centers[k].y = RandMinMax(minvy, maxvy);
                for (int kk = 0; kk < k; kk++){
                    foreach_dimension() pnt_dist.x = centers[kk].x - centers[k].x;
                    if (mynorm(pnt_dist) <= R[kk] + R[k]) {
                        flag = 0;
                        break;
                    }
                }
                fprintf(ferr, "k=%d R=%g x=%g y=%g\n", k, R[k], centers[k].x, centers[k].y);
                trials++;
                if (trials > 10000) {
                    fprintf(ferr, "ERROR: Trials are more than 10000\n");
                    break;
                }
                if (flag) {
                  fprintf(ferr, "k=%d R=%g x=%g y=%g\n", k, R[k], centers[k].x, centers[k].y);
                  break;
                }
            }
            if (i == 0) xmin_center += centers[k].x;
            if (i == Ncx - 1) xmax_center += centers[k].x;
        }
    }
    xmin_center /= Ncy - (shift_y > 0);// for satur flow true Ncy is less by 1
    xmax_center /= Ncy - (shift_y > 0);
}

double geometry(double x, double y, double z){
    coord mypoint = {x,y,z};
    coord pnt_dist;
    double phi = -HUGE;
	for (int i = 0; i < Ncx*Ncy; i++){
		foreach_dimension() pnt_dist.x = mypoint.x - centers[i].x;
       	phi = max(phi, (R[i] - mynorm(pnt_dist)));
      	// phi = -10;// no solid
	}
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

FILE *popen(const char *command, const char *mode);
int pclose(FILE *stream);

event init (t = 0) {
    R = (double *) malloc(Ncx*Ncy*sizeof(double));
    if (R == NULL) {
        fprintf(stderr, "malloc failed with R\n");
        return -1;
    }
    centers = (coord *) malloc(Ncx*Ncy*sizeof(coord));
    if (centers == NULL) {
        fprintf(stderr, "malloc failed with centers\n");
        return -1;
    }
    calc_centers(centers, R);


    if (!restore (file = "restart")) {
        int iter = 0;
        do {
            iter++;
            solid_func(fs, fs_face);
//			fraction(fs, geometry(x, y, z));
			fraction(f, bubbles(x, y, z));
			fprintf(ferr, "ITER=%d\n", iter);
        }while ((iter <=6) || ((adapt_wavelet({f, fs}, (double []){feps, fseps}, maxlevel = maxlevel, minlevel=minlevel).nf != 0) && (iter <= 15)));
        fprintf(stderr, "init refinement iter=%d\n", iter);
        foreach() {
            T[] = Tin*(1 - fs[]) + fs[]*Tcyl;
            alpha_doc[] = 0;
            u.x[] = u.y[] = 0; //u_BC*f[];// 0; //u_BC;//*(1-fs[]); // penalization will work
            un[] = u.x[];
            if (rho1 || rho2){
                rhov[] = var_hom(f[], fs[], rho1, rho2, rho3);
            }
        }
        boundary ({T, alpha_doc, u, rho});
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
		boundary((scalar *){kappa, mu});
        event("vtk_file");
    }else{
        FILE *cmd;
        char result[5000];

        cmd = popen("grep \"vtk: iter_fp\" log | awk \'{print $7}\'", "r");
        if (cmd == NULL) {
            perror("popen");
            exit(EXIT_FAILURE);
        }
        int k = 0;
        while (fgets(result, sizeof(result), cmd)) {
            printf("%s", result);
            file_timesteps[k++] = atof(result);
        }

        cmd = popen("grep \"vtk: iter_fp\" log | tail -1 | awk \'{print $4}\'", "r");
        if (cmd == NULL) {
            perror("popen");
            exit(EXIT_FAILURE);
        }
        fgets(result, sizeof(result), cmd);
        iter_fp = atoi(result) + 1;
        fprintf(ferr, "Read iter_fp+1: %d\n", iter_fp);
        pclose(cmd);
    }
}

event set_dtmax (i++) {
	RELATIVE_RES_TOLERANCE = 0.01;
    DT *= 1.05;
    double maxDT = (non_saturated>0) ? 1e-3 : 1e-2;
    DT = min(DT, maxDT);
    fprintf(ferr, "set_dtmax: tnext= %g, t=%g, DT=%g, dt=%g\n", tnext, t, DT, dt);
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

    foreach_face() uf.x[] *= 1 - fs_face.x[];
    boundary((scalar *){uf});
}
event projection(i++){
    foreach_face() uf.x[] *= 1 - fs_face.x[];
    boundary((scalar *){uf});
}



event viscous_term(i++){
    TOLERANCE = TOLERANCE_V;
//    int m_bp = max(200 - 1*i, 1);
//    double mindelta = L0/pow(2, maxlevel);
//    double nu_bp = mu1/rho1;
//    eta_s = sq(m_bp*mindelta)/nu_bp;
////    eta_s = 1e+6;
//    fprintf(ferr, "m=%d, mindelta=%g, nu=%g, eta_s=%15.12g\n", m_bp, mindelta, nu_bp, eta_s);

}

event projection(i++){
    TOLERANCE = TOLERANCE_P;
}

event end_timestep(i++){
	TOLERANCE = TOLERANCE_T;
//	relative_residual_poisson = true;
}



/*
 * xmax, ymax - maximum front tip along Ox
 * xmax_prev, ymax_prev - maximum front tip at the previous time step
 * ulocal - the local velocity at the maximum tip
 * utip - the velocity of the tip = (xmax - xmax_prev)/dt, (ymax - ymax_prev)/dt
 * volume_of_resin, volume_of_gas, volume_of_solid - the volume of resin, gas, solid
 * porosity = volume_of_gas/sq(L0) - volume of gas in the computational domain
 * porositys = volume_of_gas/(sq(L0) - volume_of_solid) - efficiency
 * vel_in_resin_x, vel_in_resin_y - the velocity of the resin in the region
 * vel_in_all_region_x, vel_in_all_region_y - the mean velocity of the resin+air+solid in the region
 * alpha_doc_avg mean value degree of cure in the resin
 * mu_avg mean value of viscosity in resin
 * T_in_resin_avg mean temperature in resin
 * T_in_fluid_avg mean temperature in fluid
 * T_in_all_region_avg mean temperature in whole region
 *
 */
double xmax_prev;
double ymax_prev;
#define region_of_averaging(x,y) ((x > xmin_center + 0.5*dist_x) && (x < xmax_center - 0.5*dist_x))
//#define region_of_averaging(x,y) (1)

struct Tip {
    coord t;
    coord ulocal;
};

scalar m[];
double time_prev = 0;
event logfile(i+=50){
    double dv_all, dv_in_resin_without_solid, dv_in_resin_solid, dv_in_resin_gas, dv_in_air_without_solid;
    double porosity, porositys;
    double xmax_tip = -1e+9, volume_all = SEPS, volume_of_fluids = SEPS, volume_of_resin = SEPS;
    double volume_of_resin_out_of_solid = SEPS, volume_of_gas = SEPS, volume_of_gas_out_of_solid = SEPS, volume_of_solid = SEPS;
    double vel_in_all_region_x = 0, vel_in_all_region_y = 0;
    double vel_in_resin_x = 0, vel_in_resin_y = 0;
    double vel_in_resin_solid_x = 0, vel_in_resin_solid_y = 0;
    double vel_in_resin_gas_x = 0, vel_in_resin_gas_y = 0;
    double gradpx_in_resin_without_solid = 0, gradpx_in_air_without_solid = 0;
    double ymax=-1e+9;
    double alpha_doc_avg = 0, mu_avg = 0, T_in_resin_avg = 0, T_in_fluid_avg = 0, T_in_all_region_avg = 0;
    int psize, pid;
    MPI_Comm_size(MPI_COMM_WORLD, &psize); // get size of procs
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    struct Tip tip[psize];

    FILE * fpv;
    FILE * fpperm;
    FILE * fppermeability;

    coord utip = {0.0,0.0,0.0}, ulocal = {0.0,0.0,0.0};
    double avg = normf_weugene(u.x, fs).avg;
    double du = change_weugene (u.x, un, fs)/(avg + SEPS); //change 1) Linf  2) un = u
    fprintf(ferr, "du= %g, du/dt: %g\n", du, du/(t - time_prev + SEPS));
    if (i == 0 ){
        xmax_prev = front_x;
        ymax_prev = 0.5*L0;
        if (pid() == 0){
            fpv = fopen("statistics_volume.txt", "wb");
            fpperm = fopen("statistics_pufmu.txt", "wb");
            fppermeability = fopen("permiability.txt", "wb");
            fprintf(fppermeability,
            "i, t, gradp, x_left, p_int[imin], x_right, p_int[imax], permeability1, permeability2, "
            "xmax_tip, xmax_prev, volume_all, volume_of_fluids, volume_of_resin, volume_of_resin_out_of_solid, "
            "volume_of_gas, volume_of_gas_out_of_solid, volume_of_solid, porosity, porositys, "
            "vel_in_all_region_x, vel_in_all_region_y, sqrt(sq(vel_in_all_region_x) + sq(vel_in_all_region_y)), "
            "vel_in_resin_x, vel_in_resin_y, sqrt(sq(vel_in_resin_x) + sq(vel_in_resin_y)), "
            "vel_in_resin_solid_x, vel_in_resin_solid_y, sqrt(sq(vel_in_resin_solid_x) + sq(vel_in_resin_solid_y)), "
            "vel_in_resin_gas_x, vel_in_resin_gas_y, sqrt(sq(vel_in_resin_gas_x) + sq(vel_in_resin_gas_y)), "
            "ulocal.x, ulocal.y, mynorm(ulocal), "
            "utip.x, utip.y, T_in_resin_avg, T_in_fluid_avg, T_in_all_region_avg, "
            "alpha_doc_avg, mu_avg, gradpx_in_resin_without_solid, gradpx_in_air_without_solid, vol_max\n");
        }
    }else{
        if (pid() == 0){
            fpv = fopen("statistics_volume.txt", "a");
            fpperm = fopen("stat_p_uf_mu.txt", "a");
            fppermeability = fopen("permiability.txt", "a");
        }
    }
    for (int j = 0; j < psize; j++){
        foreach_dimension() {
            tip[j].t.x = tip[j].ulocal.x = 0;
        }
    }
    /*
     * Calculate volumes
     */
    foreach(reduction(+:volume_all) reduction(+:volume_of_fluids) reduction(+:volume_of_resin) reduction(+:volume_of_resin_out_of_solid)
            reduction(+:volume_of_gas) reduction(+:volume_of_gas_out_of_solid) reduction(+:volume_of_solid)
            reduction(+:T_in_resin_avg) reduction(+:T_in_fluid_avg) reduction(+:T_in_all_region_avg)
            reduction(max:xmax_tip)) {
        if (f[-1] != f[1]) {
            if (x > xmax_tip) {
              xmax_tip = x;
              foreach_dimension() {
                  tip[pid].t.x = x;
                  tip[pid].ulocal.x = u.x[];
              }
//              printf("xmax_tip=%g\n", xmax_tip);
            }
        }
        if (region_of_averaging(x,y)){
            dv_all = dv();
            dv_in_resin_without_solid = f[]*(1 - fs[])*dv_all; // only in resin, exclude gas and solid !!!PHYSICAL, no penetration assumption
            dv_in_resin_solid = f[]*dv_all;      // in resin and solid, exclude gas
            dv_in_resin_gas = (1 - fs[])*dv_all; // in resin and gas, exclude solid !!! Linear velocity
            dv_in_air_without_solid = (1 - fs[])*(1 - f[])*dv_all;

            volume_all += dv_all;
            volume_of_fluids += dv_in_resin_gas;
            volume_of_resin += dv_in_resin_solid;
            volume_of_resin_out_of_solid += dv_in_resin_without_solid; // PHYSICAL, no penetration assumption
            volume_of_gas += (1 - f[])*dv_all;
            volume_of_gas_out_of_solid += (1 - fs[])*(1 - f[])*dv_all; // PHYSICAL, no penatration assumption
            volume_of_solid += fs[]*dv_all;

            T_in_resin_avg += T[]*dv_in_resin_without_solid;
            T_in_fluid_avg += T[]*dv_in_resin_gas;
            T_in_all_region_avg += T[]*dv_all;
        }
    }
    foreach(reduction(+:vel_in_all_region_x) reduction(+:vel_in_all_region_y)
            reduction(+:vel_in_resin_x) reduction(+:vel_in_resin_y)
            reduction(+:vel_in_resin_solid_x) reduction(+:vel_in_resin_solid_y)
            reduction(+:vel_in_resin_gas_x) reduction(+:vel_in_resin_gas_y)
            reduction(+:alpha_doc_avg) reduction(+:mu_avg)
            reduction(+:gradpx_in_resin_without_solid) reduction(+:gradpx_in_air_without_solid)
            ) {
        if (region_of_averaging(x,y)){
            dv_all = dv();
            dv_in_resin_without_solid = f[]*(1 - fs[])*dv_all; // only in resin, exclude gas and solid !!!PHYSICAL, no penetration assumption
            dv_in_resin_solid = f[]*dv_all;      // in resin and solid, exclude gas
            dv_in_resin_gas = (1 - fs[])*dv_all; // in resin and gas, exclude solid !!! Linear velocity
            dv_in_air_without_solid = (1 - f[])*(1 - fs[])*dv_all;

            vel_in_all_region_x += u.x[]*dv_all;
            vel_in_all_region_y += u.y[]*dv_all;
            vel_in_resin_x += u.x[]*dv_in_resin_without_solid;
            vel_in_resin_y += u.y[]*dv_in_resin_without_solid;
            vel_in_resin_solid_x += u.x[]*dv_in_resin_solid;
            vel_in_resin_solid_y += u.y[]*dv_in_resin_solid;
            vel_in_resin_gas_x += u.x[]*dv_in_resin_gas;
            vel_in_resin_gas_y += u.y[]*dv_in_resin_gas;

            alpha_doc_avg += alpha_doc[]*dv_in_resin_without_solid;
            mu_avg += 0.5*(mu.x[]+mu.y[])*dv_in_resin_without_solid;
            gradpx_in_resin_without_solid += ((p[1] - p[-1])/(2*Delta))*dv_in_resin_without_solid;
            gradpx_in_air_without_solid += ((p[1] - p[-1])/(2*Delta))*dv_in_air_without_solid;
        }
    }
    /**
    When using MPI we need to perform a global reduction to get the
    velocities and positions of the tip of a front which span multiple processes. */
    #if _MPI
        MPI_Allreduce (MPI_IN_PLACE, tip, 2*dimension*psize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    #endif

    for (int j = 0; j < psize; j++){
        if (fabs(tip[j].t.x - xmax_tip) < 1e-10) {
            ymax = tip[j].t.y;
            foreach_dimension() ulocal.x = tip[j].ulocal.x;
        }
    }

    utip.x = (i==0) ? Uin : (xmax_tip - xmax_prev)/dt; utip.y =  (i==0) ? Uin : (ymax - ymax_prev)/dt;
    vel_in_all_region_x /= volume_all; vel_in_all_region_y /= volume_all;
    vel_in_resin_x /= volume_of_resin_out_of_solid; vel_in_resin_y /= volume_of_resin_out_of_solid;
    vel_in_resin_solid_x /= volume_of_resin; vel_in_resin_solid_y /= volume_of_resin;
    vel_in_resin_gas_x /= volume_of_fluids; vel_in_resin_gas_y /= volume_of_fluids;
    porosity = volume_of_gas/volume_all;
    porositys = volume_of_gas_out_of_solid/(volume_all - volume_of_solid);
    mu_avg /= volume_of_resin_out_of_solid;
    alpha_doc_avg /= volume_of_resin_out_of_solid;
    T_in_resin_avg /= volume_of_resin_out_of_solid;
    T_in_fluid_avg /= volume_of_fluids;
    T_in_all_region_avg /= volume_all;
    gradpx_in_resin_without_solid /= volume_of_resin_out_of_solid;
    gradpx_in_air_without_solid /= volume_of_gas_out_of_solid;

    fprintf(ferr, "t= %g i= %d xmax_tip= %g xmax_prev= %g\n"
                  "volume_all= %g volume_of_fluids= %g volume_of_resin= %g volume_of_resin_out_of_solid= %g\n"
                  "volume_of_gas= %g volume_of_gas_out_of_solid= %g volume_of_solid= %g porosity= %g porositys= %g\n"
                  "vel_in_all_region_x= %g vel_in_all_region_y= %g velregion_mean= %g\n"
                  "vel_in_resin_X= %g vel_in_resin_Y= %g vel_in_resin_mean= %g\n"
                  "vel_in_resin_solid_x= %g vel_in_resin_solid_y= %g vel_in_resin_solid_mean= %g\n"
                  "vel_in_resin_gas_x= %g vel_in_resin_gas_y= %g vel_in_resin_gas_mean= %g\n"
                  "ulocal_X= %g ulocal_Y= %g ulocal_mean= %g\n"
                  "utip.x= %g utip.y= %g\n"
                  "T_in_resin_avg= %g T_in_fluid_avg= %g T_in_all_region_avg= %g\n"
                  "alpha_doc_avg= %g mu_avg= %g, gradpx_in_resin_without_solid= %g gradpx_in_air_without_solid= %g\n",
                  t, i, xmax_tip, xmax_prev,
                  volume_all, volume_of_fluids, volume_of_resin, volume_of_resin_out_of_solid,
                  volume_of_gas, volume_of_gas_out_of_solid, volume_of_solid, porosity, porositys,
                  vel_in_all_region_x, vel_in_all_region_y, sqrt(sq(vel_in_all_region_x) + sq(vel_in_all_region_y)),
                  vel_in_resin_x, vel_in_resin_y, sqrt(sq(vel_in_resin_x) + sq(vel_in_resin_y)),
                  vel_in_resin_solid_x, vel_in_resin_solid_y, sqrt(sq(vel_in_resin_solid_x) + sq(vel_in_resin_solid_y)),
                  vel_in_resin_gas_x, vel_in_resin_gas_y, sqrt(sq(vel_in_resin_gas_x) + sq(vel_in_resin_gas_y)),
                  ulocal.x, ulocal.y, mynorm(ulocal),
                  utip.x, utip.y,
                  T_in_resin_avg, T_in_fluid_avg, T_in_all_region_avg,
                  alpha_doc_avg, mu_avg, gradpx_in_resin_without_solid, gradpx_in_air_without_solid);
//    scalar m[];
    foreach() m[] = (((1 - f[])*(fs[] <=0))*region_of_averaging(x,y) > 1e-3); // m is 0 and 1 array
    boundary((scalar *){m});

    int n = tag (m); // m is modified filled be indices
    /**
    Once each cell is tagged with a unique droplet index, we can easily
    compute the volume *v* and position *b* of each droplet. Note that
    we use *foreach_leaf()* rather than *foreach()* to avoid doing a
    parallel traversal when using OpenMP. This is because we don't have
    reduction operations for the *v* and *b* arrays (yet).
     */

    double v[n];
    coord b[n];
    for (int j = 0; j < n; j++)
    v[j] = b[j].x = b[j].y = b[j].z = 0.;
    double vol_sum = 0;

    foreach_leaf()
        if (m[] > 0) {
            int j = m[] - 1;
            v[j] += dv()*m[];
            vol_sum += v[j];
            coord p = {x,y,z};
            foreach_dimension()
                b[j].x += dv()*m[]*p.x;
        }
    /**
    When using MPI we need to perform a global reduction to get the
    volumes and positions of bubbles which span multiple processes. */

    #if _MPI
        MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    #endif
    /**
    Finally we output the volume and position of each bubble to
    standard output. */
    int i_vol_max=0;
    double vol_max=-1e+9;
    for (int j = 0; j < n; j++){
        if (v[j] > vol_max) {
            vol_max = v[j];
            i_vol_max = j;
        }
        if (v[j]>0 && pid() == 0) fprintf (fpv, "i= %d t= %g j= %d v= %g bx/v= %g by/v= %g\n",
        i, t, j, v[j], b[j].x/v[j], b[j].y/v[j]);
    }
    fprintf(ferr, "i= %d t= %g i_vol_max= %d vol_max= %g\n", i, t, i_vol_max, vol_max);

    xmax_prev = xmax_tip;
    ymax_prev = ymax;

    /**
     * Permeability calculation
     * qi =
     *
     */
// Integrate to the top

    int nn = (int) pow(2, minlevel) + 1;
    double p_int[nn], uf_mean[nn][2], mu_mean[nn], xtemp, ff;

    for (int i = 0; i < nn; i++) {
        p_int[i] = 0.0;
        uf_mean[i][0] = 0.0;
        uf_mean[i][1] = 0.0;
        mu_mean[i] = 0.0;
    }

    foreach_face(x) {
        ff = 0.5*(f[] + f[-1]);
        for (int i = 0; i < nn; i++) {
            xtemp = X0 + i*L0/pow(2, minlevel);
            if (fabs(x - xtemp) < 1e-8){
              // p_int[i] += 0.5*(p[-1] + p[])*Delta*ff/L0;
              p_int[i] += p[]*Delta*f[]/L0;
              uf_mean[i][0] += uf.x[]*Delta*ff/L0;
              uf_mean[i][1] += uf.x[]*Delta*(1.0 - ff)/L0;
              mu_mean[i] += mu.x[]*ff*Delta/L0;
            }
        }
    }
    #if _MPI
        MPI_Allreduce (MPI_IN_PLACE, p_int, nn, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce (MPI_IN_PLACE, mu_mean, nn, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce (MPI_IN_PLACE, uf_mean, 2*nn, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    #endif
    double xx, dmin = 1e+9, dmax = 1e+9;
    double xmin = cyl_x + cyl_diam + 0.5*dist_x, x_left=-1e+9; // between cylinders left
    double xmax = cyl_x + cyl_diam - 0.5*dist_x + (Ncx - 1)*dist_x, x_right=1e+9; // between cylinders right
    int imin=0, imax=0;
    double dd1, dd2, gradp;
    for (int i0 = 0; i0 < nn; i0++) {
        xx = X0 + i0*L0/pow(2, minlevel);
        dd1 = fabs(xx - xmin);
        dd2 = fabs(xx - xmax);
        if (dmin > dd1){
            dmin = dd1;
            imin = i0;
            x_left = xx;
        }
        if (dmax > dd2){
            dmax = dd2;
            imax = i0;
            x_right = xx;
        }
        if (pid() == 0) fprintf(fpperm, "i0= %d t= %g x= %g p= %g mu= %g uf= %g %g\n",
                         i0, t, xx, p_int[i0], mu_mean[i0], uf_mean[i0][0], uf_mean[i0][1]);
    }
    gradp = (p_int[imax] - p_int[imin])/(dist_x*(Ncx-2));//? how to consider capillary?

    double permeability1 = mu_avg*uf_mean[imax][0]/fabs(gradpx_in_resin_without_solid + SEPS);
    double permeability2 = mu2   *uf_mean[imax][1]/fabs(gradpx_in_air_without_solid + SEPS);
    fprintf(ferr, "xmin=%g x_left=%g dmin=%g imin=%d\n", xmin, x_left, dmin, imin);
    fprintf(ferr, "xmax=%g x_right=%g dmax=%g imax=%d\n", xmax, x_right, dmax, imax);
    fprintf(ferr, "t= %g x_left= %g p_left= %g x_right= %g p_right= %g mean_dp/dx= %g permeability:K1= %g K2= %g\n",
                    t, x_left, p_int[imin], x_right, p_int[imax], gradp, permeability1, permeability2);
    if (pid() == 0){
        fprintf(fppermeability, "%d, %g, %g, %g, %g, %g, %g, %g, %g, %g, "
                            "%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, "
                            "%g, %g, %g, %g, %g, %g, %g, %g, "
                            "%g, %g, %g, %g, %g, %g, %g, %g, "
                            "%g, %g, %g, %g, %g, %g, %g, %g %g\n",
                            i, t, gradp,
                            x_left, p_int[imin],
                            x_right, p_int[imax],
                            permeability1, permeability2,
                            xmax_tip, xmax_prev,
                            volume_all, volume_of_fluids, volume_of_resin, volume_of_resin_out_of_solid, //4
                            volume_of_gas, volume_of_gas_out_of_solid, volume_of_solid, porosity, porositys, //5
                            vel_in_all_region_x, vel_in_all_region_y, sqrt(sq(vel_in_all_region_x) + sq(vel_in_all_region_y)), //3
                            vel_in_resin_x, vel_in_resin_y, sqrt(sq(vel_in_resin_x) + sq(vel_in_resin_y)), //3
                            vel_in_resin_solid_x, vel_in_resin_solid_y, sqrt(sq(vel_in_resin_solid_x) + sq(vel_in_resin_solid_y)), //3
                            vel_in_resin_gas_x, vel_in_resin_gas_y, sqrt(sq(vel_in_resin_gas_x) + sq(vel_in_resin_gas_y)), //3
                            ulocal.x, ulocal.y, mynorm(ulocal), //3
                            utip.x, utip.y, //2
                            T_in_resin_avg, T_in_fluid_avg, T_in_all_region_avg, //3
                            alpha_doc_avg, mu_avg, gradpx_in_resin_without_solid, gradpx_in_air_without_solid, vol_max); //5
        fclose(fpv);
        fclose(fpperm);
        fclose(fppermeability);
    }
    if ((non_saturated == 0) && (i > 1) && (du/(t - time_prev + SEPS) < 1e+0)) {
        fprintf (ferr, "Converged!!!\n");
        event("vtk_file");
        return 9;
    }
    time_prev = t;
    if ((p_int[imax] != p_int[imax]) || (p_int[imin] != p_int[imin]) || (vel_in_all_region_x != vel_in_all_region_x)){
        fprintf (ferr, "NAN in values!!!\n");
        event("vtk_file");
        return 9;
    }
}


// event vtk_file (i += 1)
event vtk_file (t += dt_vtk)
{
    char subname[80]; sprintf(subname, "heat_pol");
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
//    scalar curvature_by_sigma_cell[]; curvature (f, curvature_by_sigma_cell, f.sigma, add = false);
#ifdef DEBUG_BRINKMAN_PENALIZATION
    output_vtu_MPI(subname, (iter_fp) ? t + dt : 0, (scalar *) {T, alpha_doc, p, fs, f, l,  rhov}, (vector *) {u, dbp, total_rhs, a, residual_of_u, conv_term, mu, kappa});
#else
    fprintf(ferr, "output_vtu_MPI\n");
    //    output_vtu_MPI(subname, (iter_fp) ? t + dt : 0, (scalar *) {T, alpha_doc, p, fs, f}, (vector *) {u});
    output_vtu_MPI(subname, (iter_fp) ? t + dt : 0, list = (scalar *) {T, dpdx, alpha_doc, p, fs, f, l, rhov, mu_cell, m},
    vlist = (vector *) {u, a});//, mu_cell_minus, mu_cell_plus
#endif
}

//#if DUMP
event snapshot (i += snapshot_i)
//event snapshot (t += 1e-1)
//event snapshot (i += 1)
{
    char name[80];
    if (t==0)
        sprintf(name, "dump-0.0");
    else
        sprintf(name, "dump-%04g",t);
    vorticity (u, omega);
    p.nodump = false;
    dump (file = name);
}
//#endif

/**
We adapt according to the error on the embedded geometry, velocity and
tracer fields. */

#define ADAPT_SCALARS {f, fs, u.x, u.y, T, alpha_doc}
#define ADAPT_EPS_SCALARS {feps, fseps, ueps, ueps, Teps, aeps}
//#define ADAPT_SCALARS {f, fs}
//#define ADAPT_EPS_SCALARS {feps, fseps}


event adapt (i++){
	double eps_arr[] = ADAPT_EPS_SCALARS;
	MinMaxValues((scalar *) ADAPT_SCALARS, eps_arr);
	adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
//	fraction(fs, geometry(x, y, z));
    solid_func(fs, fs_face);
	//if (i > 300) stokes = true;
    foreach_face() uf.x[] *= 1 - fs_face.x[];
    boundary((scalar *){uf});
    count_cells(t, i);
}

event stop(t = 10*L0 / Uin);
