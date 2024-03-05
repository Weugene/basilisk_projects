//No centered.h
#include "grid/octree.h"
#include "run.h"
#include "timestep.h"
#include "poisson.h"
#include "utils.h"
#include "lambda2.h"
#define DEBUG_MINMAXVALUES
#include "utils-weugene.h"
#include "output_htg.h"
#include <ctype.h>
#include <dirent.h> // Directory entry
scalar fs[], f[], p[];
vector u[];
scalar l[], omega[], l2[];
double rho1=1, rho2=1, mu1=1, mu2=1;
double eta_s=1e-5; 
int NITERMIN, NITERMAX;
int iter_fp = 0;
bool relative_residual_poisson=false, relative_residual_viscous=false;
attribute {
  double sigma;
}

#ifndef rho
    #define rho(f) (clamp(f,0.,1.)*(rho1 - rho2) + rho2)
#endif
//end No centered.h


double length_min = 1e+30, length_max = -1e+30, length = 1;
double myt = 0, shiftm = 30, shiftp = 30;
char dump_name[30];

//Channel cross section Lyy*Lzz
double Vd, Vdst, deq, dst = 0.2, rst = 0.1, r_bub, l_bub;
double RhoR, MuR;
double Rho1, Rho2;
double Mu1, Mu2;
double Sigma;
double diam_tube;
double dt_vtk;
int bubcase = 22;

typedef struct {
    double Ca;    // Ca=mu*Ub/sigma
    double Vd;    // volume
    double Uc;    // Umean
    double Ud;    // Ububble
    double delta; // thickness
} Cases;
#define zv {0,0,0,0,0}
//         CA       V_d[l]    Uc[m/s]  Ud[m/s] delta*[-]
Cases cases[28]={
        zv, // 0 is empty
        zv,zv,zv,zv,zv,zv,zv, // 1-7 Air-glycerol
        {6.44e-4,  0.0349e-9, 0.0454, 0.0533, 0.1055},// 8  Air-Water delta?
        {8.42e-4,  0.0439e-9, 0.0593, 0.0697, 0.074}, // 9  Air-Water delta?
        {7.34e-4,  0.0602e-9, 0.0557, 0.0607, 0.027}, // 10 Air-Water delta?
        {6.64e-4,  0.0745e-9, 0.0534, 0.0550, 0.001}, // 11 Air-Water delta?
        {6.697e-4, 0.1893e-9, 0.0543, 0.0554, 0.006}, // 12 Air-Water
        zv,zv,zv,zv,zv, // 13-17 Air-glycerol
        {0.003, 0.1751e-9, 0.242,  0.261, 0.013},  // 18 Air-Water
        {0.008, 0.1715e-9, 0.666,  0.704, 0.023},  // 19 Air-Water
        {0.0098, 0.2208e-9, 0.757, 0.815, 0.025}, // 20 Air-Water
        {0.015,  0.1882e-9, 1.118 , 1.293, 0.039},  // 21 Air-Water
        {0.024,  0.2179e-9, 1.580, 1.944, 0.054},  // 22 Air-Water
        {0.034,  0.2179e-9, 2.060, 2.511, 1e-9},  // 23 Air-Water Ud, delta - garbage
        {0.0455, 0.2179e-9, 2.575, 3.165, 1e-9},  // 24 Air-Water Ud, delta - garbage
        {0.056,  0.2179e-9, 3.09,  4.65,   1e-9},  // 25 Air-Water Ud, delta - garbage
        {0.065,  0.2179e-9, 3.602, 5.40,   1e-9},  // 26 Air-Water Ud, delta - garbage
        {0.074,  0.2179e-9, 4.117, 6.17,   1e-9}  // 27 Air-Water Ud, delta - garbage Paper: 3
};

double Ca; // Ca = Mu*Ud/sigma
double Ca_mod; // Ca_mod = Mu*Umean/sigma
double Re; //Reynolds
double G;
double Umean;
double lDomain=30, x_init = 2;
int maxlevel_init = 0;
int maxlevel = 8;
int minlevel = 5;
int LEVEL = 6;
int adapt_method = 1; // 0 - traditional, 1 - using limitation, 2 - using array for maxlevel
int snapshot_i = 100;
double fseps = 1e-3, ueps = 1e-2;
double TOLERANCE_P = 1e-5, TOLERANCE_V = 1e-5;
double timesteps[1000];
int timesteps_count = 0;

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

// Comparison function for qsort
int compare(const void *a, const void *b) {
    double diff = *(double*)a - *(double*)b;
    return (diff < 0) ? -1 : (diff > 0);
}

void get_timesteps(){
    DIR *d; // Directory stream
    struct dirent *dir; // Pointer for directory entry

    // Open the current directory
    d = opendir(".");
    if (d) {
        while ((dir = readdir(d)) != NULL) { // Read each directory entry
            // Check if the file name starts with "dump-"
            if (strncmp(dir->d_name, "dump-", 5) == 0) {
                double time = get_double(dir->d_name);
                timesteps[timesteps_count++] = time;
            }
        }
        closedir(d); // Close the directory
        // Sort the times array
        qsort(timesteps, timesteps_count, sizeof(double), compare);

        // Print the sorted times
        for (int i = 0; i < timesteps_count; i++) {
            printf("sorted time=%f\n", timesteps[i]);
        }
    } else {
        perror("opendir() failed");
    }

}


int main (int argc, char * argv[]) {
    // set which dump files will be converted: each $(i_take)th
    // by default each dump will be converted
    if (argc > 1)
        strcpy(dump_name, argv[1]);
    else
        return 1;
    if (argc > 2)
        iter_fp = atoi (argv[2]);
    if (argc > 3)
        bubcase = atoi (argv[3]);
    if (argc > 4)
        maxlevel = atoi (argv[4]);
    DT=1e-3;
    eta_s = 1e-5;
    TOLERANCE = 1e-6;
    NITERMIN = 1;
    NITERMAX = 100;
    myt = get_double(dump_name);
    if ((bubcase >= 8 && bubcase <= 12) || bubcase >= 18) { // Water
        Rho1 = 997, Rho2 = 1.204;
        Mu1 = 0.88e-3, Mu2 = 0.019e-3;
        Sigma = 72.8e-3;
        diam_tube = 514e-6;
        dt_vtk = 1e-1;
        lDomain = (bubcase <= 21) ? 20 : 30;
    }else { // Glycerol
        Rho1 = 1250, Rho2 = 1.204;
        Mu1 = 550e-3, Mu2 = 0.019e-3;
        Sigma = 63.4e-3;
        diam_tube = 494e-6;
        dt_vtk = 1e-2;
    }
    fprintf(ferr, "dump_name=%s shiftm=%g shiftp=%g iter_fp=%d myt=%g\n", dump_name, shiftm, shiftp, iter_fp, myt);
    Ca = cases[bubcase].Ca; //Ca = Ud*Mu1/sigma
    Vd = cases[bubcase].Vd; // m^3
    Umean = cases[bubcase].Uc; // m/s
    deq = pow(6.0*Vd/pi, 1./3.);// 0.0005301091821 m
    dst = deq/diam_tube;// 1.0730955104
    rst = 0.5*dst;
    Vdst = (4./3.)*pi*cube(rst);
//    Umean = G*sq(0.5)/(8*Mu1);
    Ca_mod = Mu1*Umean/Sigma;
    Re = Umean*diam_tube*Rho1/Mu1;
    G = 32.0*Mu1*Umean/sq(diam_tube);

    r_bub = min(rst, 0.4);
    l_bub = cube(rst) / sq(r_bub);
    x_init = 1.7*l_bub;
    lDomain = 30;

    size (lDomain);
    origin (0., -L0/2., -L0/2.);
    init_grid (1 << LEVEL);

    fprintf(ferr,"BP:             eta_s=%g,     DT=%g\n"
                 "Solver:         NITERMIN=%d   NITERMAX=%d      TOLERANCE=%g  relative_residual_poisson=%d relative_residual_viscous=%d\n"
                 "OUTPUT:         dt_vtk=%g number of procs=%d\n"
                 "ADAPT:          minlevel=%d,  maxlevel=%d      adapt_meth=%d fseps=%g ueps=%g\n"
                 "Bubble case: %d\n"
                 "Properties(SI): Mu1=%g Mu2=%g Rho1=%g Rho2=%g  Sigma=%g G=%g Umean=%g\n"
                 "Apparatus:      diam_tube=%g  tube_length=%g\n"
                 "Bubble:         Vd=%g deq=%g\n",
                 eta_s, DT,
                 NITERMIN, NITERMAX, TOLERANCE, relative_residual_poisson, relative_residual_viscous,
                 dt_vtk, npe(),
                 minlevel, maxlevel, adapt_method, fseps, ueps,
                 bubcase,
                 Mu1, Mu2, Rho1, Rho2, Sigma, G, Umean,
                 diam_tube, L0,
                 Vd, deq);
    // Dimensionless parameters:
    // Averaging on diam_tube=1 and Umean=1, Mu1=1 and Rho1=1 p' = p/(Rho1*Umean^2)
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
                 "Bubble:                   Vdst=%g dst=%g  rst=%g  r_bub=%g l_bub=%g x_init=%g\n",
            mu1, mu2, rho1, rho2, f.sigma, G, Umean,
            Re, Ca, Ca_mod,
            Vdst, dst, rst, r_bub, l_bub, x_init);
    run();
}

event init (t = 0) {
    bool success = restore (file = dump_name);
    foreach(reduction(max:maxlevel_init)){
        maxlevel_init = max(maxlevel_init, level);
    }
    maxlevel = (maxlevel == 0) ? maxlevel_init : maxlevel;
    count_cells(myt, i);
    fprintf(ferr, "file has been read: L0=%g maxlevel_init=%d\n", L0, maxlevel_init);
    if (!success) {
        fprintf(ferr, "can't open the file %s. Missing this file, go to the next file\n", dump_name);
        return 0;
    }
    get_timesteps();
}

#define ADAPT_SCALARS {fs, f, p, l, omega, l2, u}
#define ADAPT_EPS_SCALARS {1e-4, 1e-4, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2}
event adapt (i++){
    if (maxlevel_init != maxlevel) {
        double eps_arr[] = ADAPT_EPS_SCALARS;
        MinMaxValues (ADAPT_SCALARS, eps_arr);
        int it = 1;
        astats s;
        do {
            fprintf(ferr, "iteration=%d\n", it);
            s = adapt_wavelet (
                    slist=ADAPT_SCALARS,
                    max=(double[]) ADAPT_EPS_SCALARS,
                    maxlevel=maxlevel,
                    minlevel=1
            );
            fprintf(ferr, "Adaptation: nf=%d nc=%d\n", s.nf, s.nc);
            if (s.nf == 0  || it > 2) break;
            it++;
        } while(1);
    }
    fprintf(ferr, "adapt with maxlevel=%d maxdepth=%d\n", maxlevel, grid->maxdepth);
}

void calculate_aux_fields(vector u, scalar l, scalar omega, scalar l2){
    foreach() l[] = level;
    vorticity (u, omega);
    lambda2 (u, l2);
}

event vtk_file (i++){
    calculate_aux_fields(u, l, omega, l2);
    vector gradp[];

    double xcg = 0, Ek = 0, Ekl = 0, Ekg = 0, volume = 0, volumeg = 0;
    length_min = 1e+30, length_max = -1e+30, length = 0;
    foreach(
        reduction(+:xcg) reduction(+:volume) reduction(+:volumeg)
        reduction(+:Ek) reduction(+:Ekl) reduction(+:Ekg)
    ) {
        if (fs[]<1){
            double dvv = (1.0 - fs[])*dv();
            double dvg = (1.0 - f[])*dvv; // gas volume
            volumeg += dvg;//gas liquid
            volume += dvv;//channel volume
            xcg   += x*dvg;// Along x
            double magu = norm(u); // speed
            Ek += 0.5*rho(f[])*dvv*sq(magu);
            Ekl += 0.5*rho1*f[]*dvv*sq(magu);
            Ekg += 0.5*rho2*(1.0 - f[])*dvv*sq(magu);
        }
        foreach_dimension() gradp.x[] = (p[] - p[-1])/Delta;
    }
    scalar sf[];
    filter_scalar_N_times(f, sf, 10);
    double volumel = volume - volumeg;
    xcg /= volumeg;
    Ek /= volume;
    Ekl /= volumel;
    Ekg /= volumeg;
    length_min = xcg - shiftm;
    length_max = xcg + shiftp;
    length = length_max - length_min;

    fprintf (
        ferr, "x= %g length_min= %g length_max= %g length= %g it_fp= %d volume= %g volumeg= %g Ek= %g Ekl= %g Ekg= %g\n",
        xcg, length_min, length_max, length, iter_fp, volume, volumeg, Ek, Ekl, Ekg
    );

    char path[] = "res"; // no slash at the end!!
    char prefix[] = "dump2pvd_compressed";
    output_htg((scalar *){fs, f, l, l2, omega, p, sf}, (vector *){u, gradp}, path, prefix, iter_fp, myt);
    fprintf(ferr, "ended output_htg\n");
    count_cells(myt, i);
    return 0;
}

event stop(t = 100); // t = 100 should  be sufficiently big in order to reach this event
