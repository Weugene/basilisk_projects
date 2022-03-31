#define BRINKMAN_PENALIZATION 1
#define DEBUG_MINMAXVALUES
#define DEBUG_BRINKMAN_PENALIZATION
#define DEBUG_MODE_TENSION
#define DEBUG_MODE_POISSON
#define DEBUG_MULTIGRID
#define REACTION_MODEL REACTION_MODEL_NON_AUTOCATALYTIC
//#define IGNORE_SOLID_MU
//#define DEBUG_OUTPUT_VTU_MPI
#define FILTERED
#define JACOBI 1
#define CORRECT_UF_FLUXES 1
#define STICKY_SOLID 1
#define CURV_PARTSTR 1// use Petr Karnakov's module
#define RELATIVE_RES
#define EPS_MAXA 2 // adapt based on Max-Min value
//#define PRINT_ALL_VALUES
#ifdef DEBUG_MODE_TENSION
scalar f_hat[];
#endif
//#define STOKES
//#define T_DIRICHLET_BC 1

scalar mu_cell[];
scalar Phi_visc[], Phi_src[];
scalar which_meth[];
scalar un[];
face vector fs_face[];
face vector av[];
static coord vel_s = {0, 0, 0};
#if CURV_PARTSTR==1
    #include "curvature_partstr.h"
#endif
#include "centered-weugene.h"
#include "rheology_model.h"
#include "tension.h"
#include "utils-weugene.h"
//#include "output_vtu_foreach.h"
#include "output_htg.h"
#include "tag.h"

int snapshot_i = 1000;
double snapshot_t = 0.5;
double dt_vtk = 0.1;
double Uin, Tin, T_solid, Tam;
double timeend = 20;
double MuRS;
coord Ggrav_ndim;
double cyl_diam, domain_size, dist_x, dist_y, cyl_x, front_x, Rbmin, Rbmax;
double dev_r, develx, devely;
double shift_x, shift_y;
int non_saturated;
bool gravityModule = false;
int Ncx, Ncy; //number of cylinders along Ox, Oy
int Nb; //number of bubbles
double layer_velocity, layer_heat;
double x_init = 2;
int maxlevel = 9;
int minlevel = 5;
int LEVEL = 7;
double maxDT;
double mindelta;
double mu_max = 0, nu_max = 0;
int adapt_method = 1; // 0 - traditional, 1 - using limitation, 2 - using array for maxlevel
double feps = 1e-10, ueps = 1e-2, rhoeps = 1e-10, Teps = 3e-2, aeps = 3e-2, mueps=1e-2;
double TOLERANCE_P = 1e-7, TOLERANCE_V = 1e-8, TOLERANCE_T = 1e-7;
double *R = NULL;
coord *centers = NULL;
char subname[150], logname[200];

#include "config_reader.h"
int main(int argc, char * argv[]) {

    struct input_yaml *input = read_config_and_assign_global_vars(argc, argv);

    MuRS = input->nums.MuRS;
    Ncy += (shift_y > 0); //for saturated flow we add artificial cylinder
    gravityModule = (fabs(Ggrav_ndim.x) + fabs(Ggrav_ndim.y) + fabs(Ggrav_ndim.z) > SEPS);
#ifdef STOKES
    stokes = true;
    stokes_heat = false;
#endif
#if T_DIRICHLET_BC != 0
    fprintf(ferr, "T_DIRICHLET_BC: 1\n");
#endif
    size(domain_size);

    const scalar temp_cyl[] = T_solid;
    T_target = temp_cyl;
    const vector U_sol[] = {vel_s.x, vel_s.y, vel_s.z};
    target_U = U_sol;
    init_grid(1 << LEVEL);
    origin (-L0/2, -L0/2.);
    periodic(top);

    a = av;
    fs.refine = fs.prolongation = fraction_refine;
    f.refine = f.prolongation = fraction_refine;

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
    partstr_conf.nohf = false; // false: both GHF and particle curvature are possible, true: use only particle
#endif
    run();
#if CURV_PARTSTR==1
    DumpCsvFin();
#endif
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
    front_x = (non_saturated) ? front_x : 1e+10;
    double limMin = X0 + Rbmax, limMax = min(min(front_x, cyl_x), X0 + L0) - Rbmax;
    double *R = malloc(max(Nb,1) * sizeof(double));
    if (R == NULL) {
        fprintf(stderr, "malloc failed with R\n");
        return -1;
    }
    coord *centers = malloc(max(Nb,1) * sizeof(coord));
    if (centers == NULL) {
        fprintf(stderr, "malloc failed with centers\n");
        return -1;
    }
// generating of Radii and centers of bubbles which are not overlapping
    srand (10); // for consistency
    int iter = 0, i = 0;
    while(i < Nb && Nb > 0){
        R[i] = RandMinMax(Rbmin, Rbmax);
        centers[i].x = RandMinMax(limMin, limMax);
        centers[i].y = RandMinMax(Y0 + Rbmax, Y0 + L0 - Rbmax);
#if dimension>2
        centers[i].z = RandMinMax(Z0 + Rbmax, Z0 + L0 - Rbmin);
#endif
        for (int j = 0; j < i; j++) {
            foreach_dimension() pnt_dist.x = centers[i].x - centers[j].x;
            if ( mynorm(pnt_dist) < 1.1*(R[i] + R[j]) ) {
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
    return min(phi, front_x - x); // with front
//    return phi; // no front
//    return (non_saturated > 0)? front_x - x : 1; // with front
}
double xmin_center = 0, xmax_center = 0;
void calc_centers(coord * centers, double * R){
    if (Ncx*Ncy == 0)
        return;
    int k;
    coord pnt_dist;
    double minvx, maxvx, minvy, maxvy;
    double limMin = cyl_x + 0.5*cyl_diam;
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
//                fprintf(ferr, "k=%d R=%g x=%g y=%g\n", k, R[k], centers[k].x, centers[k].y);
                trials++;
                if (trials > 10000) {
                    fprintf(ferr, "ERROR: Trials are more than 10000\n");
                    break;
                }
                if (flag) {
//                  fprintf(ferr, "k=%d R=%g x=%g y=%g\n", k, R[k], centers[k].x, centers[k].y);
                  break;
                }
            }
            if (i == 0) xmin_center += centers[k].x;
            if (i == Ncx - 1) xmax_center += centers[k].x;
        }
    }
    int nny = Ncy - (shift_y > 0);
    if (nny > 0)
        xmin_center /= nny;// for satur flow true Ncy is less by 1
    if (nny > 0)
        xmax_center /= nny;
}

double geometry(double x, double y, double z){
    coord mypoint = {x,y,z};
    coord pnt_dist;
    double phi = -HUGE;
	for (int i = 0; i < Ncx*Ncy; i++){
		foreach_dimension() pnt_dist.x = mypoint.x - centers[i].x;
       	phi = max(phi, (R[i] - mynorm(pnt_dist)));
//      	 phi = -10;// no solid
	}
    return phi;
}
void solid_func(scalar fs, face vector fs_face){
    vertex scalar phi[];
    foreach_vertex() {
        phi[] = geometry(x, y, z);
    }
    fractions (phi, fs, fs_face);
}

event init (t = 0) {
    char name[300];
    sprintf (name, "restart_%s", subname);
    R = (double *) malloc(max(Ncx*Ncy,1)*sizeof(double));
    if (R == NULL) {
        fprintf(stderr, "malloc failed with R\n");
        return -1;
    }
    centers = (coord *) malloc(max(Ncx*Ncy,1)*sizeof(coord));
    if (centers == NULL) {
        fprintf(stderr, "malloc failed with centers\n");
        return -1;
    }
    calc_centers(centers, R);
    if (!restore (file = name)) {
        fprintf(ferr, "The file %s can not be successfully read! Initial conditions are set\n", name);
        int it = 0;
        scalar f_smoothed[], fs_smoothed[];
        do {
            solid_func(fs, fs_face);
            fraction(f, bubbles(x, y, z));
            filter_scalar(f, f_smoothed);
            filter_scalar(fs, fs_smoothed);
        }while (adapt_wavelet({f_smoothed, fs_smoothed}, (double []){feps, feps}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
        foreach() {
            T[] = var_hom(f[], fs[], Tin, Tam, T_solid);
            alpha_doc[] = 0;
            u.x[] = u_BC*(1 - fs[]); //u_BC*f[];// 0; // penalization will work
            u.y[] = 0;
            un[] = u.x[];
            if (rho1 || rho2){
                rhov[] = var_hom(f[], fs[], rho1, rho2, rho3);
            }
        }
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
        event("vtk_file");
    }
    else
    {
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
        if (cmd == NULL) {
            fprintf(ferr, "Error in opening log file and searching iter_fp");
            perror("popen");
            exit(EXIT_FAILURE);
        }
        fgets(result, sizeof(result), cmd);
        iter_fp = atoi(result) + 1;
        fprintf(ferr, "Read last iter_fp: %d, new iter_fp: %d\n", atoi(result), iter_fp);
        pclose(cmd);
    }
}

void calc_mu(scalar muv){
    foreach(reduction(max:mu_max)){
        muv[] = mu(f[], fs[], alpha_doc[], T[]);
        if( (muv[] > mu_max) && (fs[] < 1)) mu_max = muv[];
    }
}

event properties(i++){
    solid_func(fs, fs_face);
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

event acceleration (i++) {
  if (gravityModule){
  	foreach_face()	av.x[] = Ggrav_ndim.x;
  }
}

event advection_term(i++){
    TOLERANCE = TOLERANCE_P;
}
double m_bp_iter = 30;
event viscous_term(i++){
    m_bp_iter -= 2;
    nu_max = mu_max/rho1;
    TOLERANCE = TOLERANCE_V;
    eta_s = sq(max(m_bp, m_bp_iter)*mindelta)/nu_max;
    eta_T = sq(m_bp_T*mindelta)/chi_conductivity;
   fprintf(ferr, "VISC: m=%g, mindelta=%g, nu=%g, chi_conductivity=%g, eta_s=%15.12g, eta_T=%15.12g\n", m_bp, mindelta, nu_max, chi_conductivity, eta_s, eta_T);
}

event projection(i++){
    TOLERANCE = TOLERANCE_P;
}

event end_timestep(i++){
	TOLERANCE = TOLERANCE_T;
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
//event logfile(i+=50){
//    double dv_all, dv_in_resin_without_solid, dv_in_resin_solid, dv_in_resin_gas, dv_in_air_without_solid;
//    double porosity, porositys;
//    double xmax_tip = -1e+9, volume_all = SEPS, volume_of_fluids = SEPS, volume_of_resin = SEPS;
//    double volume_of_resin_out_of_solid = SEPS, volume_of_gas = SEPS, volume_of_gas_out_of_solid = SEPS, volume_of_solid = SEPS;
//    double vel_in_all_region_x = 0, vel_in_all_region_y = 0;
//    double vel_in_resin_x = 0, vel_in_resin_y = 0;
//    double vel_in_resin_solid_x = 0, vel_in_resin_solid_y = 0;
//    double vel_in_resin_gas_x = 0, vel_in_resin_gas_y = 0;
//    double gradpx_in_resin_without_solid = 0, gradpx_in_air_without_solid = 0;
//    double ymax=-1e+9;
//    double alpha_doc_avg = 0, mu_avg = 0, T_in_resin_avg = 0, T_in_fluid_avg = 0, T_in_all_region_avg = 0;
//    int psize, pid;
//    MPI_Comm_size(MPI_COMM_WORLD, &psize); // get size of procs
//    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
//    struct Tip tip[psize];
//
//    FILE * fpv;
//    FILE * fpperm;
//    FILE * fppermeability;
//
//    coord utip = {0.0,0.0,0.0}, ulocal = {0.0,0.0,0.0};
//    double avg = normf_weugene(u.x, fs).avg;
//    double du = change_weugene (u.x, un, fs)/(avg + SEPS); //change 1) Linf  2) un = u
//    fprintf(ferr, "du= %g, du/dt: %g\n", du, du/(t - time_prev + SEPS));
//    if (i == 0 ){
//        xmax_prev = front_x;
//        ymax_prev = 0.5*L0;
//        if (pid() == 0){
//            fpv = fopen("statistics_volume.txt", "wb");
//            fpperm = fopen("statistics_pufmu.txt", "wb");
//            fppermeability = fopen("permiability.txt", "wb");
//            fprintf(fppermeability,
//            "i, t, gradp, x_left, p_int[imin], x_right, p_int[imax], permeability1, permeability2, "
//            "xmax_tip, xmax_prev, volume_all, volume_of_fluids, volume_of_resin, volume_of_resin_out_of_solid, "
//            "volume_of_gas, volume_of_gas_out_of_solid, volume_of_solid, porosity, porositys, "
//            "vel_in_all_region_x, vel_in_all_region_y, sqrt(sq(vel_in_all_region_x) + sq(vel_in_all_region_y)), "
//            "vel_in_resin_x, vel_in_resin_y, sqrt(sq(vel_in_resin_x) + sq(vel_in_resin_y)), "
//            "vel_in_resin_solid_x, vel_in_resin_solid_y, sqrt(sq(vel_in_resin_solid_x) + sq(vel_in_resin_solid_y)), "
//            "vel_in_resin_gas_x, vel_in_resin_gas_y, sqrt(sq(vel_in_resin_gas_x) + sq(vel_in_resin_gas_y)), "
//            "ulocal.x, ulocal.y, mynorm(ulocal), "
//            "utip.x, utip.y, T_in_resin_avg, T_in_fluid_avg, T_in_all_region_avg, "
//            "alpha_doc_avg, mu_avg, gradpx_in_resin_without_solid, gradpx_in_air_without_solid, vol_max\n");
//        }
//    }else{
//        if (pid() == 0){
//            fpv = fopen("statistics_volume.txt", "a");
//            fpperm = fopen("stat_p_uf_mu.txt", "a");
//            fppermeability = fopen("permiability.txt", "a");
//        }
//    }
//    for (int j = 0; j < psize; j++){
//        foreach_dimension() {
//            tip[j].t.x = tip[j].ulocal.x = 0;
//        }
//    }
//    /*
//     * Calculate volumes
//     */
//    foreach(reduction(+:volume_all) reduction(+:volume_of_fluids) reduction(+:volume_of_resin) reduction(+:volume_of_resin_out_of_solid)
//            reduction(+:volume_of_gas) reduction(+:volume_of_gas_out_of_solid) reduction(+:volume_of_solid)
//            reduction(+:T_in_resin_avg) reduction(+:T_in_fluid_avg) reduction(+:T_in_all_region_avg)
//            reduction(max:xmax_tip)) {
//        if (f[-1] != f[1]) {
//            if (x > xmax_tip) {
//              xmax_tip = x;
//              foreach_dimension() {
//                  tip[pid].t.x = x;
//                  tip[pid].ulocal.x = u.x[];
//              }
//    //              printf("xmax_tip=%g\n", xmax_tip);
//            }
//        }
//        if (region_of_averaging(x,y)){
//            dv_all = dv();
//            dv_in_resin_without_solid = f[]*(1 - fs[])*dv_all; // only in resin, exclude gas and solid !!!PHYSICAL, no penetration assumption
//            dv_in_resin_solid = f[]*dv_all;      // in resin and solid, exclude gas
//            dv_in_resin_gas = (1 - fs[])*dv_all; // in resin and gas, exclude solid !!! Linear velocity
//            dv_in_air_without_solid = (1 - fs[])*(1 - f[])*dv_all;
//
//            volume_all += dv_all;
//            volume_of_fluids += dv_in_resin_gas;
//            volume_of_resin += dv_in_resin_solid;
//            volume_of_resin_out_of_solid += dv_in_resin_without_solid; // PHYSICAL, no penetration assumption
//            volume_of_gas += (1 - f[])*dv_all;
//            volume_of_gas_out_of_solid += (1 - fs[])*(1 - f[])*dv_all; // PHYSICAL, no penatration assumption
//            volume_of_solid += fs[]*dv_all;
//
//            T_in_resin_avg += T[]*dv_in_resin_without_solid;
//            T_in_fluid_avg += T[]*dv_in_resin_gas;
//            T_in_all_region_avg += T[]*dv_all;
//        }
//    }
//    foreach(reduction(+:vel_in_all_region_x) reduction(+:vel_in_all_region_y)
//            reduction(+:vel_in_resin_x) reduction(+:vel_in_resin_y)
//            reduction(+:vel_in_resin_solid_x) reduction(+:vel_in_resin_solid_y)
//            reduction(+:vel_in_resin_gas_x) reduction(+:vel_in_resin_gas_y)
//            reduction(+:alpha_doc_avg) reduction(+:mu_avg)
//            reduction(+:gradpx_in_resin_without_solid) reduction(+:gradpx_in_air_without_solid)
//            ) {
//        if (region_of_averaging(x,y)){
//            dv_all = dv();
//            dv_in_resin_without_solid = f[]*(1 - fs[])*dv_all; // only in resin, exclude gas and solid !!!PHYSICAL, no penetration assumption
//            dv_in_resin_solid = f[]*dv_all;      // in resin and solid, exclude gas
//            dv_in_resin_gas = (1 - fs[])*dv_all; // in resin and gas, exclude solid !!! Linear velocity
//            dv_in_air_without_solid = (1 - f[])*(1 - fs[])*dv_all;
//
//            vel_in_all_region_x += u.x[]*dv_all;
//            vel_in_all_region_y += u.y[]*dv_all;
//            vel_in_resin_x += u.x[]*dv_in_resin_without_solid;
//            vel_in_resin_y += u.y[]*dv_in_resin_without_solid;
//            vel_in_resin_solid_x += u.x[]*dv_in_resin_solid;
//            vel_in_resin_solid_y += u.y[]*dv_in_resin_solid;
//            vel_in_resin_gas_x += u.x[]*dv_in_resin_gas;
//            vel_in_resin_gas_y += u.y[]*dv_in_resin_gas;
//
//            alpha_doc_avg += alpha_doc[]*dv_in_resin_without_solid;
//            mu_avg += 0.5*(mu.x[]+mu.y[])*dv_in_resin_without_solid;
//            gradpx_in_resin_without_solid += ((p[1] - p[-1])/(2*Delta))*dv_in_resin_without_solid;
//            gradpx_in_air_without_solid += ((p[1] - p[-1])/(2*Delta))*dv_in_air_without_solid;
//        }
//    }
//    /**
//    When using MPI we need to perform a global reduction to get the
//    velocities and positions of the tip of a front which span multiple processes. */
//    #if _MPI
//        MPI_Allreduce (MPI_IN_PLACE, tip, 2*dimension*psize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//    #endif
//
//    for (int j = 0; j < psize; j++){
//        if (fabs(tip[j].t.x - xmax_tip) < 1e-10) {
//            ymax = tip[j].t.y;
//            foreach_dimension() ulocal.x = tip[j].ulocal.x;
//        }
//    }
//    utip.x = (i==0) ? Uin : (xmax_tip - xmax_prev)/dt; utip.y =  (i==0) ? Uin : (ymax - ymax_prev)/dt;
//    vel_in_all_region_x /= volume_all; vel_in_all_region_y /= volume_all;
//    vel_in_resin_x /= volume_of_resin_out_of_solid; vel_in_resin_y /= volume_of_resin_out_of_solid;
//    vel_in_resin_solid_x /= volume_of_resin; vel_in_resin_solid_y /= volume_of_resin;
//    vel_in_resin_gas_x /= volume_of_fluids; vel_in_resin_gas_y /= volume_of_fluids;
//    porosity = volume_of_gas/volume_all;
//    porositys = volume_of_gas_out_of_solid/(volume_all - volume_of_solid);
//    mu_avg /= volume_of_resin_out_of_solid;
//    alpha_doc_avg /= volume_of_resin_out_of_solid;
//    T_in_resin_avg /= volume_of_resin_out_of_solid;
//    T_in_fluid_avg /= volume_of_fluids;
//    T_in_all_region_avg /= volume_all;
//    gradpx_in_resin_without_solid /= volume_of_resin_out_of_solid;
//    gradpx_in_air_without_solid /= volume_of_gas_out_of_solid;
//
//    fprintf(ferr, "t= %g i= %d xmax_tip= %g xmax_prev= %g\n"
//                  "volume_all= %g volume_of_fluids= %g volume_of_resin= %g volume_of_resin_out_of_solid= %g\n"
//                  "volume_of_gas= %g volume_of_gas_out_of_solid= %g volume_of_solid= %g porosity= %g porositys= %g\n"
//                  "vel_in_all_region_x= %g vel_in_all_region_y= %g velregion_mean= %g\n"
//                  "vel_in_resin_X= %g vel_in_resin_Y= %g vel_in_resin_mean= %g\n"
//                  "vel_in_resin_solid_x= %g vel_in_resin_solid_y= %g vel_in_resin_solid_mean= %g\n"
//                  "vel_in_resin_gas_x= %g vel_in_resin_gas_y= %g vel_in_resin_gas_mean= %g\n"
//                  "ulocal_X= %g ulocal_Y= %g ulocal_mean= %g\n"
//                  "utip.x= %g utip.y= %g\n"
//                  "T_in_resin_avg= %g T_in_fluid_avg= %g T_in_all_region_avg= %g\n"
//                  "alpha_doc_avg= %g mu_avg= %g, gradpx_in_resin_without_solid= %g gradpx_in_air_without_solid= %g\n",
//                  t, i, xmax_tip, xmax_prev,
//                  volume_all, volume_of_fluids, volume_of_resin, volume_of_resin_out_of_solid,
//                  volume_of_gas, volume_of_gas_out_of_solid, volume_of_solid, porosity, porositys,
//                  vel_in_all_region_x, vel_in_all_region_y, sqrt(sq(vel_in_all_region_x) + sq(vel_in_all_region_y)),
//                  vel_in_resin_x, vel_in_resin_y, sqrt(sq(vel_in_resin_x) + sq(vel_in_resin_y)),
//                  vel_in_resin_solid_x, vel_in_resin_solid_y, sqrt(sq(vel_in_resin_solid_x) + sq(vel_in_resin_solid_y)),
//                  vel_in_resin_gas_x, vel_in_resin_gas_y, sqrt(sq(vel_in_resin_gas_x) + sq(vel_in_resin_gas_y)),
//                  ulocal.x, ulocal.y, mynorm(ulocal),
//                  utip.x, utip.y,
//                  T_in_resin_avg, T_in_fluid_avg, T_in_all_region_avg,
//                  alpha_doc_avg, mu_avg, gradpx_in_resin_without_solid, gradpx_in_air_without_solid);
//    //    scalar m[];
//    foreach() m[] = (((1 - f[])*(fs[] <=0))*region_of_averaging(x,y) > 1e-3); // m is 0 and 1 array
//
//    int n = tag (m); // m is modified filled be indices
//    /**
//    Once each cell is tagged with a unique droplet index, we can easily
//    compute the volume *v* and position *b* of each droplet. Note that
//    we use *foreach_leaf()* rather than *foreach()* to avoid doing a
//    parallel traversal when using OpenMP. This is because we don't have
//    reduction operations for the *v* and *b* arrays (yet).
//     */
//
//    double v[n];
//    coord b[n];
//    for (int j = 0; j < n; j++)
//    v[j] = b[j].x = b[j].y = b[j].z = 0.;
//    double vol_sum = 0;
//
//        foreach_leaf()
//        if (m[] > 0) {
//            int j = m[] - 1;
//            v[j] += dv()*m[];
//            vol_sum += v[j];
//            coord p = {x,y,z};
//            foreach_dimension()
//                b[j].x += dv()*m[]*p.x;
//        }
//    /**
//    When using MPI we need to perform a global reduction to get the
//    volumes and positions of bubbles which span multiple processes. */
//
//    #if _MPI
//        MPI_Allreduce (MPI_IN_PLACE, v, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//        MPI_Allreduce (MPI_IN_PLACE, b, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//    #endif
//    /**
//    Finally we output the volume and position of each bubble to
//    standard output. */
//    int i_vol_max=0;
//    double vol_max=-1e+9;
//    for (int j = 0; j < n; j++){
//        if (v[j] > vol_max) {
//            vol_max = v[j];
//            i_vol_max = j;
//        }
//        if (v[j]>0 && pid() == 0) fprintf (fpv, "i= %d t= %g j= %d v= %g bx/v= %g by/v= %g\n",
//        i, t, j, v[j], b[j].x/v[j], b[j].y/v[j]);
//    }
//    fprintf(ferr, "i= %d t= %g i_vol_max= %d vol_max= %g\n", i, t, i_vol_max, vol_max);
//
//    xmax_prev = xmax_tip;
//    ymax_prev = ymax;
//
//    /**
//     * Permeability calculation
//     * qi =
//     *
//     */
//    // Integrate to the top
//
//    int nn = (int) pow(2, minlevel) + 1;
//    double p_int[nn], uf_mean[nn][2], mu_mean[nn], xtemp, ff;
//
//    for (int i = 0; i < nn; i++) {
//        p_int[i] = 0.0;
//        uf_mean[i][0] = 0.0;
//        uf_mean[i][1] = 0.0;
//        mu_mean[i] = 0.0;
//    }
//
//    foreach_face(x) {
//        ff = 0.5*(f[] + f[-1]);
//        for (int i = 0; i < nn; i++) {
//            xtemp = X0 + i*L0/pow(2, minlevel);
//            if (fabs(x - xtemp) < 1e-8){
//              // p_int[i] += 0.5*(p[-1] + p[])*Delta*ff/L0;
//              p_int[i] += p[]*Delta*f[]/L0;
//              uf_mean[i][0] += uf.x[]*Delta*ff/L0;
//              uf_mean[i][1] += uf.x[]*Delta*(1.0 - ff)/L0;
//              mu_mean[i] += mu.x[]*ff*Delta/L0;
//            }
//        }
//    }
//    #if _MPI
//        MPI_Allreduce (MPI_IN_PLACE, p_int, nn, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//        MPI_Allreduce (MPI_IN_PLACE, mu_mean, nn, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//        MPI_Allreduce (MPI_IN_PLACE, uf_mean, 2*nn, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
//    #endif
//    double xx, dmin = 1e+9, dmax = 1e+9;
//    double xmin = cyl_x + cyl_diam + 0.5*dist_x, x_left=-1e+9; // between cylinders left
//    double xmax = cyl_x + cyl_diam - 0.5*dist_x + (Ncx - 1)*dist_x, x_right=1e+9; // between cylinders right
//    int imin=0, imax=0;
//    double dd1, dd2, gradp;
//    for (int i0 = 0; i0 < nn; i0++) {
//        xx = X0 + i0*L0/pow(2, minlevel);
//        dd1 = fabs(xx - xmin);
//        dd2 = fabs(xx - xmax);
//        if (dmin > dd1){
//            dmin = dd1;
//            imin = i0;
//            x_left = xx;
//        }
//        if (dmax > dd2){
//            dmax = dd2;
//            imax = i0;
//            x_right = xx;
//        }
//        if (pid() == 0) fprintf(fpperm, "i0= %d t= %g x= %g p= %g mu= %g uf= %g %g\n",
//                         i0, t, xx, p_int[i0], mu_mean[i0], uf_mean[i0][0], uf_mean[i0][1]);
//    }
//    gradp = (p_int[imax] - p_int[imin])/(dist_x*(Ncx-2));//? how to consider capillary?
//
//    double permeability1 = mu_avg*uf_mean[imax][0]/fabs(gradpx_in_resin_without_solid + SEPS);
//    double permeability2 = mu2   *uf_mean[imax][1]/fabs(gradpx_in_air_without_solid + SEPS);
//    fprintf(ferr, "xmin=%g x_left=%g dmin=%g imin=%d\n", xmin, x_left, dmin, imin);
//    fprintf(ferr, "xmax=%g x_right=%g dmax=%g imax=%d\n", xmax, x_right, dmax, imax);
//    fprintf(ferr, "t= %g x_left= %g p_left= %g x_right= %g p_right= %g mean_dp/dx= %g permeability:K1= %g K2= %g\n",
//                    t, x_left, p_int[imin], x_right, p_int[imax], gradp, permeability1, permeability2);
//    if (pid() == 0){
//        fprintf(fppermeability, "%d, %g, %g, %g, %g, %g, %g, %g, %g, %g, "
//                            "%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, "
//                            "%g, %g, %g, %g, %g, %g, %g, %g, "
//                            "%g, %g, %g, %g, %g, %g, %g, %g, "
//                            "%g, %g, %g, %g, %g, %g, %g, %g %g\n",
//                            i, t, gradp,
//                            x_left, p_int[imin],
//                            x_right, p_int[imax],
//                            permeability1, permeability2,
//                            xmax_tip, xmax_prev,
//                            volume_all, volume_of_fluids, volume_of_resin, volume_of_resin_out_of_solid, //4
//                            volume_of_gas, volume_of_gas_out_of_solid, volume_of_solid, porosity, porositys, //5
//                            vel_in_all_region_x, vel_in_all_region_y, sqrt(sq(vel_in_all_region_x) + sq(vel_in_all_region_y)), //3
//                            vel_in_resin_x, vel_in_resin_y, sqrt(sq(vel_in_resin_x) + sq(vel_in_resin_y)), //3
//                            vel_in_resin_solid_x, vel_in_resin_solid_y, sqrt(sq(vel_in_resin_solid_x) + sq(vel_in_resin_solid_y)), //3
//                            vel_in_resin_gas_x, vel_in_resin_gas_y, sqrt(sq(vel_in_resin_gas_x) + sq(vel_in_resin_gas_y)), //3
//                            ulocal.x, ulocal.y, mynorm(ulocal), //3
//                            utip.x, utip.y, //2
//                            T_in_resin_avg, T_in_fluid_avg, T_in_all_region_avg, //3
//                            alpha_doc_avg, mu_avg, gradpx_in_resin_without_solid, gradpx_in_air_without_solid, vol_max); //5
//        fclose(fpv);
//        fclose(fpperm);
//        fclose(fppermeability);
//    }
//    if ((non_saturated == 0) && (i > 1) && (du/(t - time_prev + SEPS) < 1e+0)) {
//        fprintf (ferr, "Converged!!!\n");
//        event("vtk_file");
//        return 9;
//    }
//    time_prev = t;
//    if ((p_int[imax] != p_int[imax]) || (p_int[imin] != p_int[imin]) || (vel_in_all_region_x != vel_in_all_region_x)){
//        fprintf (ferr, "NAN in values!!!\n");
//        event("vtk_file");
//        return 9;
//    }
//}


//event vtk_file (i += 1)
//event vtk_file (t += dt_vtk)
//event vtk_file (i += 1000)
//{
//    char name[300];
//    sprintf (name, "vtk_%s", subname);
//    scalar l[], dpdx[];
//    foreach() {
//        l[] = level;
//        dpdx[] = (p[1] - p[-1])/(2*Delta);
//    }
//    calcPhiVisc (u, uf, T, alpha_doc, f, fs, Phi_visc, Phi_src); //TODO: visc dissipation?
//
//#ifdef DEBUG_BRINKMAN_PENALIZATION
//    output_vtu_MPI(name, (iter_fp) ? t + dt : 0, list = (scalar *) {T, alpha_doc, p, dpdx, fs, f, l, rhov, mu_cell, my_kappa, which_meth, Phi_visc, Phi_src},
//    vlist = (vector *) {u, g, uf, av, dbp, total_rhs, residual_of_u, divtauu, fs_face});
//#else
//    output_vtu_MPI(name, (iter_fp) ? t + dt : 0, list = (scalar *) {T, dpdx, alpha_doc, p, fs, f, l, rhov, mu_cell, m, Phi_visc, Phi_src},
//    vlist = (vector *) {u, av});
//#endif
//}


event report(i+=1000){
    char path[]="res"; // no slash at the end!!
    char prefix[] = "data";
    scalar l[], dpdx[];
    foreach() {
        l[] = level;
        dpdx[] = (p[1] - p[-1])/(2*Delta);
    }
    calcPhiVisc (u, uf, T, alpha_doc, f, fs, Phi_visc, Phi_src); //TODO: visc dissipation?

    output_htg(path, prefix, (iter_fp) ? t + dt : 0, (scalar *) {T, alpha_doc, p, dpdx, fs, f, l, rhov, mu_cell, my_kappa, which_meth, Phi_visc, Phi_src},
           (vector *){u, g, uf, av, dbp, total_rhs, residual_of_u, divtauu, fs_face});
}

/**
We adapt according to the error on the embedded geometry, velocity and
tracer fields. */

#define ADAPT_SCALARS {rhov, fs, u.x, u.y, T, alpha_doc, mu_cell}
#define ADAPT_EPS_SCALARS {rhoeps, feps, ueps, ueps, Teps, aeps, mueps}

event adapt (i++){
	double eps_arr[] = ADAPT_EPS_SCALARS;
	MinMaxValues((scalar *) ADAPT_SCALARS, eps_arr);
	adapt_wavelet ((scalar *) ADAPT_SCALARS, eps_arr, maxlevel = maxlevel, minlevel = minlevel);
    if (i%100==0) count_cells(t, i);
}

//event snapshot (i += snapshot_i)
event snapshot (t += snapshot_t)
{
    char name[300];
    if (t==0)
        sprintf (name, "dump_%s-0.0", subname);
    else
        sprintf (name, "dump_%s-%g", subname, t);
    dump (file = name);
}

event check_fail(i += 100){
    foreach(serial, noauto){
        if ((u.x[] != u.x[]) || (fabs(u.x[]) > 10e+10)) {
            printf("Nan values: x=%g y=%g |u|=%g, pid()=%d\n", x, y, fabs(u.x[]), pid());
            return 9;
        }
    }
}

event stop(t = timeend);
