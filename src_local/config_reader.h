/******************************************************************************
* CYAML schema to tell libcyaml about both expected YAML and data structure.
*
* (Our CYAML schema is just a bunch of static const data.)
******************************************************************************/

#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cyaml/cyaml.h>

//colors
#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

// nondimensional numbers
struct numbers {
    double Re; // Reynolds number Re = rho*U*d/mu
    double Ca; // Capillary number Ca = Mu*U*d/sigma
    double Pe; // Peclet number Pe = Cp*Rho*U*d/Kappa
    double Pr; // Prandtl number Pr = Mu*Cp/Kappa
    double Fr; //sqrt(sq(U)/(grav*d));
    double RhoR;
    double RhoRS;
    double MuR;
    double MuRS;
    double CpR;
    double CpRS;
    double KappaR;
    double KappaRS;
};

struct nondimensionalalized_vars {
    double domain_size;
//    double characteristic_size;
//    double Uin;
//    double Tin;
    double T_solid;
    double Tam;
    double sigma_ndim;
    double Ggrav_ndim[3];
    double Htr;
    double Arrhenius_const;
    double Ea_by_R;
    double n_degree;
    double Eeta_by_Rg;
    double chi;
    double rho[3]; // density
    double kappa[3]; // heat conduction
    double Cp[3]; // heat capacity
    double mu[4]; // dynamic viscosity
    double nu[4]; // kinematic viscosity
    double ratio_front_x;
    double ratio_dist_x;
    double ratio_dist_y;
    double ratio_Rbmin;
    double ratio_Rbmax;
    double cyl_x;
    double shift_x;
    double shift_y;
    double dev_r;
    double develx;
    double devely;
    int Nb;
    int Ncx;
    int Ncy;
    int non_saturated;
};

// dimensional variables
struct dimensional_vars {
    double domain_size;
    double characteristic_size;
    double Uin;
    double Tin;
    double Tam;
    double T_solid;
    double Sigma;
    double Grav[3];
    double Htr;
    double Arrhenius_const;
    double Ea_by_R;
    double n_degree;
    double Eeta_by_Rg;
    double chi;
    double Rho[3]; // density
    double Kappa[3]; // heat conduction
    double CP[3]; // heat capacity
    double Mu[4]; // heat capacity
};

struct numerical_params {
    int minlevel;
    int LEVEL;
    int maxlevel;
    double time;
    int iter_fp;
    int snapshot_i;
    double snapshot_t;
    double dt_vtk;
    int N_smooth;
    double TOLERANCE;
    double TOLERANCE_V;
    double TOLERANCE_P;
    double TOLERANCE_T;
    int NITERMIN;
    int NITERMAX;
    double CFL;
    double CFL_ARR;
    double DT;
    double maxDT;
    double m_bp;
    double m_bp_T;
    double feps;
    double ueps;
    double rhoeps;
    double Teps;
    double aeps;
    double mueps;
};

struct input_yaml {
    char* name;
    struct numbers nums;
    struct dimensional_vars dv;
    struct nondimensionalalized_vars ndv;
    struct numerical_params num_params;
};



static const cyaml_schema_value_t bar_schema = {
        CYAML_VALUE_FLOAT(CYAML_FLAG_DEFAULT, double)
};

static const cyaml_schema_field_t numbers_fields[] = {
        CYAML_FIELD_FLOAT("Re", CYAML_FLAG_OPTIONAL, struct numbers, Re),
        CYAML_FIELD_FLOAT("Ca", CYAML_FLAG_OPTIONAL, struct numbers, Ca),
        CYAML_FIELD_FLOAT("Pe", CYAML_FLAG_OPTIONAL, struct numbers, Pe),
        CYAML_FIELD_FLOAT("Pr", CYAML_FLAG_OPTIONAL, struct numbers, Pr),
        CYAML_FIELD_FLOAT("Fr", CYAML_FLAG_OPTIONAL, struct numbers, Fr),
        CYAML_FIELD_FLOAT("RhoR", CYAML_FLAG_OPTIONAL, struct numbers, RhoR),
        CYAML_FIELD_FLOAT("RhoRS", CYAML_FLAG_OPTIONAL, struct numbers, RhoRS),
        CYAML_FIELD_FLOAT("MuR", CYAML_FLAG_OPTIONAL, struct numbers, MuR),
        CYAML_FIELD_FLOAT("MuRS", CYAML_FLAG_OPTIONAL, struct numbers, MuRS),
        CYAML_FIELD_FLOAT("CpR", CYAML_FLAG_OPTIONAL, struct numbers, CpR),
        CYAML_FIELD_FLOAT("CpRS", CYAML_FLAG_OPTIONAL, struct numbers, CpRS),
        CYAML_FIELD_FLOAT("KappaR", CYAML_FLAG_OPTIONAL, struct numbers, KappaR),
        CYAML_FIELD_FLOAT("KappaRS", CYAML_FLAG_OPTIONAL, struct numbers, KappaRS),
        CYAML_FIELD_END
};

static const cyaml_schema_field_t dimensional_fields[] = {
        CYAML_FIELD_FLOAT("domain_size", CYAML_FLAG_OPTIONAL, struct dimensional_vars, domain_size),
        CYAML_FIELD_FLOAT("characteristic_size", CYAML_FLAG_DEFAULT, struct dimensional_vars, characteristic_size),
        CYAML_FIELD_FLOAT("Uin", CYAML_FLAG_DEFAULT, struct dimensional_vars, Uin),
        CYAML_FIELD_FLOAT("Tin", CYAML_FLAG_DEFAULT, struct dimensional_vars, Tin),
        CYAML_FIELD_FLOAT("Tam", CYAML_FLAG_DEFAULT, struct dimensional_vars, Tam),
        CYAML_FIELD_FLOAT("T_solid", CYAML_FLAG_DEFAULT, struct dimensional_vars, T_solid),
        CYAML_FIELD_FLOAT("Sigma", CYAML_FLAG_DEFAULT, struct dimensional_vars, Sigma),
        CYAML_FIELD_SEQUENCE_FIXED("Grav", CYAML_FLAG_DEFAULT, struct dimensional_vars, Grav, &bar_schema, 3),
        CYAML_FIELD_FLOAT("Htr", CYAML_FLAG_DEFAULT, struct dimensional_vars, Htr),
        CYAML_FIELD_FLOAT("Arrhenius_const", CYAML_FLAG_DEFAULT, struct dimensional_vars, Arrhenius_const),
        CYAML_FIELD_FLOAT("Ea_by_R", CYAML_FLAG_DEFAULT, struct dimensional_vars, Ea_by_R),
        CYAML_FIELD_FLOAT("n_degree", CYAML_FLAG_DEFAULT, struct dimensional_vars, n_degree),
        CYAML_FIELD_FLOAT("Eeta_by_Rg", CYAML_FLAG_DEFAULT, struct dimensional_vars, Eeta_by_Rg),
        CYAML_FIELD_FLOAT("chi", CYAML_FLAG_DEFAULT, struct dimensional_vars, chi),
        CYAML_FIELD_SEQUENCE_FIXED("Rho", CYAML_FLAG_DEFAULT, struct dimensional_vars, Rho, &bar_schema, 3),
        CYAML_FIELD_SEQUENCE_FIXED("Kappa", CYAML_FLAG_DEFAULT, struct dimensional_vars, Kappa, &bar_schema, 3),
        CYAML_FIELD_SEQUENCE_FIXED("CP", CYAML_FLAG_DEFAULT, struct dimensional_vars, CP, &bar_schema, 3),
        CYAML_FIELD_SEQUENCE_FIXED("Mu", CYAML_FLAG_OPTIONAL, struct dimensional_vars, Mu, &bar_schema, 4),
        CYAML_FIELD_END
};

static const cyaml_schema_field_t non_dimensional_fields[] = {
        CYAML_FIELD_FLOAT("domain_size", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, domain_size),
        CYAML_FIELD_FLOAT("T_solid", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, T_solid),
        CYAML_FIELD_FLOAT("Tam", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, Tam),
        CYAML_FIELD_FLOAT("sigma_ndim", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, sigma_ndim),
        CYAML_FIELD_SEQUENCE_FIXED("Ggrav_ndim", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, Ggrav_ndim, &bar_schema, 3),
        CYAML_FIELD_FLOAT("Htr", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, Htr),
        CYAML_FIELD_FLOAT("Arrhenius_const", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, Arrhenius_const),
        CYAML_FIELD_FLOAT("Ea_by_R", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, Ea_by_R),
        CYAML_FIELD_FLOAT("n_degree", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, n_degree),
        CYAML_FIELD_FLOAT("Eeta_by_Rg", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, Eeta_by_Rg),
        CYAML_FIELD_FLOAT("chi", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, chi),
        CYAML_FIELD_SEQUENCE_FIXED("rho", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, rho, &bar_schema, 3),
        CYAML_FIELD_SEQUENCE_FIXED("kappa", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, kappa, &bar_schema, 3),
        CYAML_FIELD_SEQUENCE_FIXED("Cp", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, Cp, &bar_schema, 3),
        CYAML_FIELD_SEQUENCE_FIXED("mu", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, mu, &bar_schema, 4),
        CYAML_FIELD_SEQUENCE_FIXED("nu", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, nu, &bar_schema, 4),
        CYAML_FIELD_FLOAT("ratio_front_x", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, ratio_front_x),
        CYAML_FIELD_FLOAT("ratio_dist_x", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, ratio_dist_x),
        CYAML_FIELD_FLOAT("ratio_dist_y", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, ratio_dist_y),
        CYAML_FIELD_FLOAT("ratio_Rbmin", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, ratio_Rbmin),
        CYAML_FIELD_FLOAT("ratio_Rbmax", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, ratio_Rbmax),
        CYAML_FIELD_FLOAT("cyl_x", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, cyl_x),
        CYAML_FIELD_FLOAT("shift_x", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, shift_x),
        CYAML_FIELD_FLOAT("shift_y", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, shift_y),
        CYAML_FIELD_FLOAT("dev_r", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, dev_r),
        CYAML_FIELD_FLOAT("develx", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, develx),
        CYAML_FIELD_FLOAT("devely", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, devely),
        CYAML_FIELD_INT("Nb", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, Nb),
        CYAML_FIELD_INT("Ncx", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, Ncx),
        CYAML_FIELD_INT("Ncy", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, Ncy),
        CYAML_FIELD_INT("non_saturated", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, non_saturated),
        CYAML_FIELD_END
};

static const cyaml_schema_field_t numerical_params_fields[] = {
        CYAML_FIELD_INT("minlevel", CYAML_FLAG_DEFAULT, struct numerical_params, minlevel),
        CYAML_FIELD_INT("LEVEL", CYAML_FLAG_DEFAULT, struct numerical_params, LEVEL),
        CYAML_FIELD_INT("maxlevel", CYAML_FLAG_DEFAULT, struct numerical_params, maxlevel),
        CYAML_FIELD_FLOAT("time", CYAML_FLAG_DEFAULT, struct numerical_params, time),
        CYAML_FIELD_INT("iter_fp", CYAML_FLAG_DEFAULT, struct numerical_params, iter_fp),
        CYAML_FIELD_INT("snapshot_i", CYAML_FLAG_OPTIONAL, struct numerical_params, snapshot_i),
        CYAML_FIELD_FLOAT("snapshot_t", CYAML_FLAG_OPTIONAL, struct numerical_params, snapshot_t),
        CYAML_FIELD_FLOAT("dt_vtk", CYAML_FLAG_OPTIONAL, struct numerical_params, dt_vtk),
        CYAML_FIELD_INT("N_smooth", CYAML_FLAG_OPTIONAL, struct numerical_params, N_smooth),
        CYAML_FIELD_FLOAT("TOLERANCE", CYAML_FLAG_OPTIONAL, struct numerical_params, TOLERANCE),
        CYAML_FIELD_FLOAT("TOLERANCE_V", CYAML_FLAG_OPTIONAL, struct numerical_params, TOLERANCE_V),
        CYAML_FIELD_FLOAT("TOLERANCE_P", CYAML_FLAG_OPTIONAL, struct numerical_params, TOLERANCE_P),
        CYAML_FIELD_FLOAT("TOLERANCE_T", CYAML_FLAG_OPTIONAL, struct numerical_params, TOLERANCE_T),
        CYAML_FIELD_INT("NITERMIN", CYAML_FLAG_OPTIONAL, struct numerical_params, NITERMIN),
        CYAML_FIELD_INT("NITERMAX", CYAML_FLAG_OPTIONAL, struct numerical_params, NITERMAX),
        CYAML_FIELD_FLOAT("CFL", CYAML_FLAG_OPTIONAL, struct numerical_params, CFL),
        CYAML_FIELD_FLOAT("CFL_ARR", CYAML_FLAG_OPTIONAL, struct numerical_params, CFL_ARR),
        CYAML_FIELD_FLOAT("maxDT", CYAML_FLAG_OPTIONAL, struct numerical_params, maxDT),
        CYAML_FIELD_FLOAT("DT", CYAML_FLAG_DEFAULT, struct numerical_params, DT),
        CYAML_FIELD_FLOAT("m_bp", CYAML_FLAG_OPTIONAL, struct numerical_params, m_bp),
        CYAML_FIELD_FLOAT("m_bp_T", CYAML_FLAG_OPTIONAL, struct numerical_params, m_bp_T),
        CYAML_FIELD_FLOAT("feps", CYAML_FLAG_OPTIONAL, struct numerical_params, feps),
        CYAML_FIELD_FLOAT("ueps", CYAML_FLAG_OPTIONAL, struct numerical_params, ueps),
        CYAML_FIELD_FLOAT("rhoeps", CYAML_FLAG_OPTIONAL, struct numerical_params, rhoeps),
        CYAML_FIELD_FLOAT("Teps", CYAML_FLAG_OPTIONAL, struct numerical_params, Teps),
        CYAML_FIELD_FLOAT("aeps", CYAML_FLAG_OPTIONAL, struct numerical_params, aeps),
        CYAML_FIELD_FLOAT("mueps", CYAML_FLAG_OPTIONAL, struct numerical_params, mueps),
        CYAML_FIELD_END
};

/* CYAML mapping schema fields array for the top level mapping. */
static const cyaml_schema_field_t top_mapping_schema[] = {
        CYAML_FIELD_STRING_PTR("name", CYAML_FLAG_POINTER,struct input_yaml, name,0, CYAML_UNLIMITED),
        CYAML_FIELD_MAPPING("dimensional_vars", CYAML_FLAG_DEFAULT,struct input_yaml, dv, dimensional_fields),
        CYAML_FIELD_MAPPING("non_dimensional_vars", CYAML_FLAG_DEFAULT,struct input_yaml, ndv, non_dimensional_fields),
        CYAML_FIELD_MAPPING("nums", CYAML_FLAG_OPTIONAL,struct input_yaml, nums, numbers_fields),
        CYAML_FIELD_MAPPING("num_params", CYAML_FLAG_DEFAULT,struct input_yaml, num_params, numerical_params_fields),
        CYAML_FIELD_END
};

/* CYAML value schema for the top level mapping. */
static const cyaml_schema_value_t top_schema = {
        CYAML_VALUE_MAPPING(CYAML_FLAG_POINTER,
        struct input_yaml, top_mapping_schema)
};


/******************************************************************************
 * Actual code to load and save YAML doc using libcyaml.
 ******************************************************************************/

/* Our CYAML config.
 *
 * If you want to change it between calls, don't make it const.
 *
 */
static const cyaml_config_t config = {
        .log_fn = cyaml_log,            /* Use the default logging function. */
        .mem_fn = cyaml_mem,            /* Use the default memory allocator. */
        .log_level = CYAML_LOG_WARNING, /* Logging errors and warnings only. CYAML_LOG_WARNING, CYAML_LOG_DEBUG */
};

struct input_yaml* read_config(int argc, char *argv[])
{
    struct input_yaml *input = (struct input_yaml *) calloc(1, sizeof(struct input_yaml));
    cyaml_err_t err;
    enum {
        ARG_PROG_NAME,
        ARG_PATH_IN,
        ARG__COUNT,
    };

    /* Handle args */
    if (argc != ARG__COUNT) {
        fprintf(stderr, "Usage:\n");
        fprintf(stderr, "  %s <INPUT>\n", argv[ARG_PROG_NAME]);
        return input;
    }

    /* Load input file. */
    err = cyaml_load_file(argv[ARG_PATH_IN], &config,
                          &top_schema, (void **)&input, NULL);
    if (err != CYAML_OK) {
        fprintf(stderr, "ERROR: %s\n", cyaml_strerror(err));
        return input;
    }


    /* Use the data. */
    fprintf(fout, "%s:\n", input->name);
    struct dimensional_vars *dv = &input->dv;
    struct nondimensionalalized_vars *ndv = &input->ndv;
    struct numbers *nums = &input->nums;
    struct numerical_params *numpar = &input->num_params;

//    Mu0*exp(Eeta_by_Rg/Tin)
    if (!input->dv.Mu[0]) dv->Mu[0] = 3.85e-7;
    if (!input->dv.Mu[1]) dv->Mu[1] = dv->Mu[0]*exp(dv->Eeta_by_Rg/dv->Tin);
    if (!input->dv.Mu[2]) dv->Mu[2] = 1.963e-2;//1.963e-5
    if (!input->dv.Mu[3]) dv->Mu[3] = dv->Mu[0]*exp(dv->Eeta_by_Rg/dv->T_solid);
    if (!input->dv.domain_size && !input->ndv.domain_size){
        dv->domain_size = dv->characteristic_size*ndv->ratio_dist_x*max( max(ndv->Ncx, ndv->Ncy),1);
        ndv->domain_size = ndv->ratio_dist_x*max( max(ndv->Ncx, ndv->Ncy),1);
    }
//    compute non dimensional nums:
    if (!input->nums.Re) nums->Re = dv->Uin*dv->characteristic_size*dv->Rho[0]/dv->Mu[1];
    if (!input->nums.Ca) nums->Ca = dv->Uin*dv->Mu[1]/dv->Sigma;
    if (!input->nums.Pe) nums->Pe = dv->Uin*dv->characteristic_size*dv->CP[0]*dv->Rho[0]/dv->Kappa[0];
    if (!input->nums.Pr) nums->Pr = dv->Mu[0]*dv->CP[0]/dv->Kappa[0];
    if (!input->nums.Fr) nums->Fr = sqrt(dv->Uin*dv->Uin/(sqrt(dv->Grav[0]*dv->Grav[0] + dv->Grav[1]*dv->Grav[1] + dv->Grav[2]*dv->Grav[2])*dv->characteristic_size + 1e-10));
    if (!input->nums.RhoR) nums->RhoR = dv->Rho[1]/dv->Rho[0];
    if (!input->nums.RhoRS) nums->RhoRS = dv->Rho[2]/dv->Rho[0];
    if (!input->nums.MuR) nums->MuR = dv->Mu[2]/dv->Mu[1];
    if (!input->nums.MuRS) nums->MuRS = dv->Mu[3]/dv->Mu[1];
    if (!input->nums.CpR) nums->CpR = dv->CP[1]/dv->CP[0];
    if (!input->nums.CpRS) nums->CpRS = dv->CP[2]/dv->CP[0];
    if (!input->nums.KappaR) nums->KappaR = dv->Kappa[1]/dv->Kappa[0];
    if (!input->nums.KappaRS) nums->KappaRS = dv->Kappa[2]/dv->Kappa[0];

    if (!input->ndv.domain_size) ndv->domain_size = dv->domain_size/dv->characteristic_size;


    ndv->T_solid = dv->T_solid/dv->Tin;
    ndv->Tam = dv->Tam/dv->Tin;
    ndv->sigma_ndim = 1.0/(nums->Re*nums->Ca + 1e-10);
    ndv->Ggrav_ndim[0] = dv->Grav[0]/(nums->Fr*nums->Fr);
    ndv->Ggrav_ndim[1] = dv->Grav[1]/(nums->Fr*nums->Fr);
    ndv->Ggrav_ndim[2] = dv->Grav[2]/(nums->Fr*nums->Fr);
    ndv->Htr = dv->Htr/(dv->CP[0]*dv->Tin);
    ndv->Arrhenius_const = dv->Arrhenius_const*dv->characteristic_size/(dv->Uin + 1e-10);
    ndv->Ea_by_R = dv->Ea_by_R/dv->Tin;
    ndv->n_degree = dv->n_degree;
    ndv->Eeta_by_Rg = dv->Eeta_by_Rg/dv->Tin;
    ndv->chi = dv->chi;
    ndv->rho[0] = 1;
    ndv->rho[1] = nums->RhoR;
    ndv->rho[2] = nums->RhoRS;
    ndv->kappa[0] = dv->Kappa[0]/(dv->Rho[0]*dv->CP[0]*dv->characteristic_size*dv->Uin + 1e-10);
    ndv->kappa[1] = ndv->kappa[0]*nums->KappaR;
    ndv->kappa[2] = ndv->kappa[0]*nums->KappaRS;
    ndv->Cp[0] = 1.0;
    ndv->Cp[1] = ndv->Cp[0]*nums->CpR;
    ndv->Cp[2] = ndv->Cp[0]*nums->CpRS;
    ndv->mu[0] = (1.0/nums->Re)*(dv->Mu[0]/dv->Mu[1]);
    ndv->mu[1] = (1.0/nums->Re);
    ndv->mu[2] = (1.0/nums->Re)*nums->MuR;
    ndv->mu[3] = (1.0/nums->Re)*nums->MuRS;
    ndv->nu[0] = ndv->mu[0]/ndv->rho[0];
    ndv->nu[1] = ndv->mu[1]/ndv->rho[0];
    ndv->nu[2] = ndv->mu[2]/ndv->rho[1];
    ndv->nu[3] = ndv->mu[3]/ndv->rho[2];

//    setting up default values for non-dimensional vars
    if (!input->ndv.ratio_front_x) ndv->ratio_front_x = 100;
    if (!input->ndv.ratio_dist_x) ndv->ratio_dist_x = 2;
    if (!input->ndv.ratio_dist_y) ndv->ratio_dist_y = 2;
    if (!input->ndv.ratio_Rbmin) ndv->ratio_Rbmin = 0.2;
    if (!input->ndv.ratio_Rbmax) ndv->ratio_Rbmax = 1;
    if (!input->ndv.cyl_x) ndv->cyl_x = 0;
    if (!input->ndv.shift_x) ndv->shift_x = 0;
    if (!input->ndv.shift_y) ndv->shift_y = 0;
    if (!input->ndv.dev_r) ndv->dev_r = 0;
    if (!input->ndv.develx) ndv->develx = 0;
    if (!input->ndv.devely) ndv->devely = 0;
    if (!input->ndv.Nb) ndv->Nb = 0;
    if (!input->ndv.Ncx) ndv->Ncx = 0;
    if (!input->ndv.Ncy) ndv->Ncy = 0;
    if (!input->ndv.non_saturated) ndv->non_saturated = 0;

//    setting up default values for numerical simulation params
    if (!input->num_params.N_smooth) numpar->N_smooth = 1;
    if (!input->num_params.snapshot_i) numpar->snapshot_i = 1000;
    if (!input->num_params.snapshot_t) numpar->snapshot_t = 0.5;
    if (!input->num_params.dt_vtk) numpar->dt_vtk = 0.1;
    if (!input->num_params.TOLERANCE) numpar->TOLERANCE = 1e-7;
    if (!input->num_params.TOLERANCE_V) numpar->TOLERANCE_V = 1e-8;
    if (!input->num_params.TOLERANCE_P) numpar->TOLERANCE_P = 1e-7;
    if (!input->num_params.TOLERANCE_T) numpar->TOLERANCE_T = 1e-7;
    if (!input->num_params.NITERMIN) numpar->NITERMIN = 1;
    if (!input->num_params.NITERMAX) numpar->NITERMAX = 100;
    if (!input->num_params.CFL) numpar->CFL = 0.5;
    if (!input->num_params.CFL_ARR) numpar->CFL_ARR = 0.5;
    if (!input->num_params.maxDT) numpar->maxDT = 1e+10;
    if (!input->num_params.m_bp) numpar->m_bp = 2.0;
    if (!input->num_params.m_bp_T) numpar->m_bp_T = 2.0;
    if (!input->num_params.feps) numpar->feps = 1e-10;
    if (!input->num_params.ueps) numpar->ueps = 1e-2;
    if (!input->num_params.rhoeps) numpar->rhoeps = 1e-10;
    if (!input->num_params.Teps) numpar->Teps = 1e-2;
    if (!input->num_params.aeps) numpar->aeps = 1e-2;
    if (!input->num_params.mueps) numpar->mueps = 1e-2;

    fprintf(fout, "%s============ Dimensional params ============\n%s", KRED, KNRM);
    fprintf(fout, "characteristic_size=%g [m] ", dv->characteristic_size);
    fprintf(fout, "Uin=%g [m/s]\n", dv->Uin);
    fprintf(fout, "Tin=%g [K] ", dv->Tin);
    fprintf(fout, "Tam=%g [K] ", dv->Tam);
    fprintf(fout, "T_solid=%g [K]\n", dv->T_solid);
    fprintf(fout, "Sigma=%g [N/m]\n", dv->Sigma);
    fprintf(fout, "Grav=[%g %g %g] [m/s^2]\n", dv->Grav[0], dv->Grav[1], dv->Grav[2]);
    fprintf(fout, "Htr=%g [J/kg] ", dv->Htr);
    fprintf(fout, "Arrhenius_const=%g [1/s] ", dv->Arrhenius_const);
    fprintf(fout, "Ea_by_R=%g [K] ", dv->Ea_by_R);
    fprintf(fout, "n_degree=%g\n", dv->n_degree);
    fprintf(fout, "Eeta_by_Rg=%g [K] ", dv->Eeta_by_Rg);
    fprintf(fout, "chi=%g\n", dv->chi);
    fprintf(fout, "Rho=[%g %g %g] [kg/m^3]\n", dv->Rho[0], dv->Rho[1], dv->Rho[2]);
    fprintf(fout, "Kappa=[%g %g %g] [W/(m*K)]\n", dv->Kappa[0], dv->Kappa[1], dv->Kappa[2]);
    fprintf(fout, "CP=[%g %g %g] [J/(kg*K)]\n", dv->CP[0], dv->CP[1], dv->CP[2]);
    fprintf(fout, "Mu=[%g %g %g %g] [Pa*s]\n", dv->Mu[0], dv->Mu[1], dv->Mu[2], dv->Mu[3]);

    fprintf(fout, "%s============ Non - dimensional numbers ============\n%s", KRED, KNRM);
    fprintf(fout, "Re=%g\n", nums->Re);
    fprintf(fout, "Ca=%g\n", nums->Ca);
    fprintf(fout, "Pe=%g\n", nums->Pe);
    fprintf(fout, "Pr=%g\n", nums->Pr);
    fprintf(fout, "Fr=%g\n", nums->Fr);
    fprintf(fout, "RhoR=%g\n", nums->RhoR);
    fprintf(fout, "RhoRS=%g\n", nums->RhoRS);
    fprintf(fout, "MuR=%g\n", nums->MuR);
    fprintf(fout, "MuRS=%g\n", nums->MuRS);
    fprintf(fout, "CpR=%g\n", nums->CpR);
    fprintf(fout, "CpRS=%g\n", nums->CpRS);
    fprintf(fout, "KappaR=%g\n", nums->KappaR);
    fprintf(fout, "KappaRS=%g\n", nums->KappaRS);

    fprintf(fout, "%s============ Non - dimensional params ============\n%s", KRED, KNRM);
    fprintf(fout, "domain_size=%g ", ndv->domain_size);
    fprintf(fout, "Tam=%g ", ndv->Tam);
    fprintf(fout, "T_solid=%g\n", ndv->T_solid);
    fprintf(fout, "sigma_ndim=%g\n", ndv->sigma_ndim);
    fprintf(fout, "Grav_ndim=[%g %g %g]\n", ndv->Ggrav_ndim[0], ndv->Ggrav_ndim[1], ndv->Ggrav_ndim[2]);
    fprintf(fout, "Htr=%g ", ndv->Htr);
    fprintf(fout, "Arrhenius_const=%g ", ndv->Arrhenius_const);
    fprintf(fout, "Ea_by_R=%g ", ndv->Ea_by_R);
    fprintf(fout, "n_degree=%g\n", ndv->n_degree);
    fprintf(fout, "Eeta_by_Rg=%g ", ndv->Eeta_by_Rg);
    fprintf(fout, "chi=%g\n", ndv->chi);
    fprintf(fout, "rho=[%g %g %g]\n", ndv->rho[0], ndv->rho[1], ndv->rho[2]);
    fprintf(fout, "kappa=[%g %g %g]\n", ndv->kappa[0], ndv->kappa[1], ndv->kappa[2]);
    fprintf(fout, "Cp=[%g %g %g]\n", ndv->Cp[0], ndv->Cp[1], ndv->Cp[2]);
    fprintf(fout, "mu=[%g %g %g %g]\n", ndv->mu[0], ndv->mu[1], ndv->mu[2], ndv->mu[3]);
    fprintf(fout, "nu=[%g %g %g %g]\n", ndv->nu[0], ndv->nu[1], ndv->nu[2], ndv->nu[3]);
    fprintf(fout, "ratio_front_x=%g ratio_dist_x=%g ratio_dist_y=%g\n", ndv->ratio_front_x, ndv->ratio_dist_x, ndv->ratio_dist_y);
    fprintf(fout, "ratio_Rbmin=%g ratio_Rbmax=%g\n", ndv->ratio_Rbmin, ndv->ratio_Rbmax);
    fprintf(fout, "cyl_x=%g shift_x=%g shift_y=%g\n", ndv->cyl_x, ndv->shift_x, ndv->shift_y);
    fprintf(fout, "dev_r=%g develx=%g devely=%g\n", ndv->dev_r, ndv->develx, ndv->devely);
    fprintf(fout, "Nb=%d Ncx=%d Ncy=%d non_saturated=%d\n", ndv->Nb, ndv->Ncx, ndv->Ncy, ndv->non_saturated);

    fprintf(fout, "%s============ Numerical simulation params ============\n%s", KRED, KNRM);
    fprintf(fout, "minlevel=%d LEVEL=%d maxlevel=%d\n", numpar->minlevel, numpar->LEVEL, numpar->maxlevel);
    fprintf(fout, "time=%g\n", numpar->time);
    fprintf(fout, "iter_fp=%d\n", numpar->iter_fp);
    fprintf(fout, "snapshot_i=%d snapshot_t=%g dt_vtk=%g\n", numpar->snapshot_i, numpar->snapshot_t, numpar->dt_vtk);
    fprintf(fout, "N_smooth=%d\n", numpar->N_smooth);
    fprintf(fout, "TOLERANCE=%g TOLERANCE_V=%g TOLERANCE_P=%g TOLERANCE_T=%g\n", numpar->TOLERANCE, numpar->TOLERANCE_V, numpar->TOLERANCE_P, numpar->TOLERANCE_T);
    fprintf(fout, "NITERMIN=%d NITERMAX=%d\n", numpar->NITERMIN, numpar->NITERMAX);
    fprintf(fout, "CFL=%g CFL_ARR=%g\n", numpar->CFL, numpar->CFL_ARR);
    fprintf(fout, "DT=%g maxDT=%g\n", numpar->DT, numpar->maxDT);
    fprintf(fout, "m_bp=%g m_bp_T=%g\n", numpar->m_bp, numpar->m_bp_T);
    fprintf(fout, "layer_velocity=%g layer_heat=%g\n", 1.0/sqrt(input->nums.Re), 1.0/sqrt(input->nums.Pe));
    fprintf(fout, "feps=%g ueps=%g rhoeps=%g Teps=%g aeps=%g mueps=%g\n", numpar->feps, numpar->ueps, numpar->rhoeps, numpar->Teps, numpar->aeps, numpar->mueps);

    /* Free the data */
//    cyaml_free(&config, &top_schema, input, 0);
    return input;
}


struct input_yaml* read_config_and_assign_global_vars(int argc, char *argv[])
{
    struct input_yaml *input = read_config(argc, argv);

    // numerical simulation params
    sprintf(subname, "%s", input->name);
    sprintf(logname, "log_%s", input->name);
    snapshot_i = input->num_params.snapshot_i;
    snapshot_t = input->num_params.snapshot_t;
    dt_vtk = input->num_params.dt_vtk;
    minlevel = input->num_params.minlevel;
    LEVEL = input->num_params.LEVEL;
    timeend = input->num_params.time;
    maxlevel = input->num_params.maxlevel;
    feps = input->num_params.feps;
    ueps = input->num_params.ueps;
    rhoeps = input->num_params.rhoeps;
    Teps = input->num_params.Teps;
    aeps = input->num_params.aeps;
    mueps = input->num_params.mueps;
    TOLERANCE = input->num_params.TOLERANCE;
    TOLERANCE_P = input->num_params.TOLERANCE_P;
    TOLERANCE_V = input->num_params.TOLERANCE_V;
    TOLERANCE_T = input->num_params.TOLERANCE_T;
    NITERMIN = input->num_params.NITERMIN;
    NITERMAX = input->num_params.NITERMAX;
    maxDT = input->num_params.maxDT;
    DT = input->num_params.DT;

    CFL = input->num_params.CFL;
    CFL_ARR = input->num_params.CFL_ARR;
    N_smooth = input->num_params.N_smooth;

    // dimensionless parameters:
    cyl_diam = 1;
    Uin = 1;
    Tin = 1;
    T_solid = input->ndv.T_solid;
    Tam = input->ndv.Tam;
    Ggrav_ndim.x = input->ndv.Ggrav_ndim[0];
    Ggrav_ndim.y = input->ndv.Ggrav_ndim[1];
    Ggrav_ndim.z = input->ndv.Ggrav_ndim[2];
    f.sigma = (fabs(input->ndv.sigma_ndim) < 1e+10) ? input->ndv.sigma_ndim : 0;
    domain_size = input->ndv.domain_size;
    dist_x = input->ndv.ratio_dist_x;
    dist_y = input->ndv.ratio_dist_y;
    front_x = input->ndv.ratio_front_x;
    Rbmin = input->ndv.ratio_Rbmin;
    Rbmax = input->ndv.ratio_Rbmax;
    cyl_x = input->ndv.cyl_x;
    shift_x = input->ndv.shift_x;
    shift_y = input->ndv.shift_y;
    dev_r = input->ndv.dev_r;
    develx = input->ndv.develx;
    devely = input->ndv.devely;
    non_saturated = input->ndv.non_saturated;
    Ncx = input->ndv.Ncx;
    Ncy = input->ndv.Ncy;
    Nb = input->ndv.Nb;
    layer_velocity = 1.0/sqrt(input->nums.Re);
    layer_heat = 1.0/sqrt(input->nums.Pe);
    m_bp = input->num_params.m_bp;
    m_bp_T = input->num_params.m_bp_T;
    domain_size = input->ndv.domain_size;
    mindelta = domain_size/pow(2, maxlevel);

    rho1 = input->ndv.rho[0];
    rho2 = input->ndv.rho[1];
    rho3 = input->ndv.rho[2];
    mu0 = input->ndv.mu[0];
    mu1 = input->ndv.mu[1];
    mu2 = input->ndv.mu[2];
    mu3 = input->ndv.mu[3];
    Cp1 = input->ndv.Cp[0];
    Cp2 = input->ndv.Cp[1];
    Cp3 = input->ndv.Cp[2];
    kappa1 = input->ndv.kappa[0];
    kappa2 = input->ndv.kappa[1];
    kappa3 = input->ndv.kappa[2];


    chi_conductivity = kappa1 / (rho1 * Cp1);
    eta_s = sq(m_bp*mindelta)/(mu1/rho1);
    eta_T = sq(m_bp_T*mindelta)/chi_conductivity;

    Htr = input->ndv.Htr;
    Arrhenius_const = input->ndv.Arrhenius_const;
    Ea_by_R = input->ndv.Ea_by_R;
    n_degree = input->ndv.n_degree;
    Eeta_by_Rg = input->ndv.Eeta_by_Rg;
    chi = input->ndv.chi;
    Arrhenius_const = input->ndv.Arrhenius_const;
    Arrhenius_const = input->ndv.Arrhenius_const;
    Arrhenius_const = input->ndv.Arrhenius_const;

    return input;
}

