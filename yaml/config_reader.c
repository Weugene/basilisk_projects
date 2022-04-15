/******************************************************************************
* CYAML schema to tell libcyaml about both expected YAML and data structure.
*
* (Our CYAML schema is just a bunch of static const data.)
******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include "config_reader.h"
#include "math.h"

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
    float Re; // Reynolds number Re = rho*U*d/mu
    float Ca; // Capillary number Ca = Mu*U*d/sigma
    float Pe; // Peclet number Pe = Cp*Rho*U*d/Kappa
    float Pr; // Prandtl number Pr = Mu*Cp/Kappa
    float Fr; //sqrt(sq(U)/(grav*d));
    float RhoR;
    float RhoRS;
    float MuR;
    float MuRS;
    float CpR;
    float CpRS;
    float KappaR;
    float KappaRS;
};

struct nondimensionalalized_vars {
    float domain_size;
//    float characteristic_size;
//    float Uin;
//    float Tin;
    float T_solid;
    float Tam;
    float sigma_ndim;
    float Ggrav_ndim[3];
    float Htr;
    float Arrhenius_const;
    float Ea_by_R;
    float n_degree;
    float Eeta_by_Rg;
    float chi;
    float rho[3]; // density
    float kappa[3]; // heat conduction
    float Cp[3]; // heat capacity
    float mu[4]; // heat capacity
    float ratio_front_x;
    float ratio_dist_x;
    float ratio_dist_y;
    float ratio_Rbmin;
    float ratio_Rbmax;
    float shift_x;
    float shift_y;
    float dev_r;
    float develx;
    float devely;
    int Nb;
    int Ncx;
    int Ncy;
    int non_saturated;
};

// dimensional variables
struct dimensional_vars {
    float domain_size;
    float characteristic_size;
    float Uin;
    float Tin;
    float Tam;
    float T_solid;
    float Sigma;
    float Grav[3];
    float Htr;
    float Arrhenius_const;
    float Ea_by_R;
    float n_degree;
    float Eeta_by_Rg;
    float chi;
    float Rho[3]; // density
    float Kappa[3]; // heat conduction
    float CP[3]; // heat capacity
    float Mu[4]; // heat capacity
};

struct numerical_params {
    int minlevel;
    int maxlevel;
    int iter_fp;
    int snapshot_i;
    float snapshot_t;
    float dt_vtk;
    int N_smooth;
    float TOLERANCE;
    float TOLERANCEV;
    float TOLERANCEP;
    int NITERMIN;
    int NITERMAX;
    float CFL;
    float CFL_ARR;
    float DT;
    float maxDT;
    float m_bp;
    float m_bp_T;
    float feps;
    float ueps;
    float rhoeps;
    float Teps;
    float aeps;
    float mueps;
};

struct input_yaml {
    char* name;
    struct numbers nums;
    struct dimensional_vars dv;
    struct nondimensionalalized_vars ndv;
    struct numerical_params num_params;
};

static const cyaml_schema_value_t bar_schema = {
        CYAML_VALUE_FLOAT(CYAML_FLAG_DEFAULT, float)
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
        CYAML_FIELD_FLOAT("ratio_front_x", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, ratio_front_x),
        CYAML_FIELD_FLOAT("ratio_dist_x", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, ratio_dist_x),
        CYAML_FIELD_FLOAT("ratio_dist_y", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, ratio_dist_y),
        CYAML_FIELD_FLOAT("ratio_Rbmin", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, ratio_Rbmin),
        CYAML_FIELD_FLOAT("ratio_Rbmax", CYAML_FLAG_OPTIONAL, struct nondimensionalalized_vars, ratio_Rbmax),
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
        CYAML_FIELD_INT("maxlevel", CYAML_FLAG_DEFAULT, struct numerical_params, maxlevel),
        CYAML_FIELD_INT("iter_fp", CYAML_FLAG_DEFAULT, struct numerical_params, iter_fp),
        CYAML_FIELD_INT("snapshot_i", CYAML_FLAG_OPTIONAL, struct numerical_params, snapshot_i),
        CYAML_FIELD_FLOAT("snapshot_t", CYAML_FLAG_OPTIONAL, struct numerical_params, snapshot_t),
        CYAML_FIELD_FLOAT("dt_vtk", CYAML_FLAG_OPTIONAL, struct numerical_params, dt_vtk),
        CYAML_FIELD_INT("N_smooth", CYAML_FLAG_OPTIONAL, struct numerical_params, N_smooth),
        CYAML_FIELD_FLOAT("TOLERANCE", CYAML_FLAG_OPTIONAL, struct numerical_params, TOLERANCE),
        CYAML_FIELD_FLOAT("TOLERANCEV", CYAML_FLAG_OPTIONAL, struct numerical_params, TOLERANCEV),
        CYAML_FIELD_FLOAT("TOLERANCEP", CYAML_FLAG_OPTIONAL, struct numerical_params, TOLERANCEP),
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
        CYAML_FIELD_MAPPING("nums", CYAML_FLAG_DEFAULT,struct input_yaml, nums, numbers_fields),
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

int read_config(int argc, char *argv[])
{
    cyaml_err_t err;
    struct input_yaml *input;
    enum {
        ARG_PROG_NAME,
        ARG_PATH_IN,
        ARG__COUNT,
    };

    /* Handle args */
    if (argc != ARG__COUNT) {
        fprintf(stderr, "Usage:\n");
        fprintf(stderr, "  %s <INPUT>\n", argv[ARG_PROG_NAME]);
        return EXIT_FAILURE;
    }

    /* Load input file. */
    err = cyaml_load_file(argv[ARG_PATH_IN], &config,
                          &top_schema, (void **)&input, NULL);
    if (err != CYAML_OK) {
        fprintf(stderr, "ERROR: %s\n", cyaml_strerror(err));
        return EXIT_FAILURE;
    }


    /* Use the data. */
    printf("%s:\n", input->name);
    struct dimensional_vars *dv = &input->dv;
    struct nondimensionalalized_vars *ndv = &input->ndv;
    struct numbers *nums = &input->nums;
    struct numerical_params *numpar = &input->num_params;

//    Mu0*exp(Eeta_by_Rg/Tin)
    if (!input->dv.Mu[0]) dv->Mu[0] = 3.85e-7;
    if (!input->dv.Mu[1]) dv->Mu[1] = dv->Mu[0]*exp(dv->Eeta_by_Rg/dv->Tin);
    if (!input->dv.Mu[2]) dv->Mu[2] = 1.963e-5;
    if (!input->dv.Mu[3]) dv->Mu[3] = dv->Mu[0]*exp(dv->Eeta_by_Rg/dv->T_solid)*dv->Rho[2]/dv->Rho[0];
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
    ndv->kappa[1] = dv->Kappa[0]*nums->KappaR;
    ndv->kappa[2] = dv->Kappa[0]*nums->KappaRS;
    ndv->Cp[0] = 1.0;
    ndv->Cp[1] = ndv->Cp[0]*nums->CpR;
    ndv->Cp[2] = ndv->Cp[0]*nums->CpRS;
    ndv->mu[0] = (1.0/nums->Re)*(dv->Mu[0]/dv->Mu[1]);
    ndv->mu[1] = (1.0/nums->Re);
    ndv->mu[2] = (1.0/nums->Re)*nums->MuR;
    ndv->mu[3] = (1.0/nums->Re)*nums->MuRS;

//    setting up default values for non-dimensional vars
    if (!input->ndv.ratio_front_x) ndv->ratio_front_x = 100;
    if (!input->ndv.ratio_dist_x) ndv->ratio_dist_x = 2;
    if (!input->ndv.ratio_dist_y) ndv->ratio_dist_y = 2;
    if (!input->ndv.ratio_Rbmin) ndv->ratio_Rbmin = 0.2;
    if (!input->ndv.ratio_Rbmax) ndv->ratio_Rbmax = 1;
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
    if (!input->num_params.TOLERANCEV) numpar->TOLERANCEV = 1e-8;
    if (!input->num_params.TOLERANCEP) numpar->TOLERANCEP = 1e-7;
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

    printf("%s============ Dimensional params ============\n%s", KRED, KNRM);
    printf("characteristic_size=%g [m] ", dv->characteristic_size);
    printf("Uin=%g [m/s]\n", dv->Uin);
    printf("Tin=%g [K] ", dv->Tin);
    printf("T_solid=%g [K]\n", dv->T_solid);
    printf("Sigma=%g [N/m]\n", dv->Sigma);
    printf("Grav=[%g %g %g] [m/s^2]\n", dv->Grav[0], dv->Grav[1], dv->Grav[2]);
    printf("Htr=%g [J/kg] ", dv->Htr);
    printf("Arrhenius_const=%g [1/s] ", dv->Arrhenius_const);
    printf("Ea_by_R=%g [K] ", dv->Ea_by_R);
    printf("n_degree=%g\n", dv->n_degree);
    printf("Eeta_by_Rg=%g [K] ", dv->Eeta_by_Rg);
    printf("chi=%g\n", dv->chi);
    printf("Rho=[%g %g %g] [kg/m^3]\n", dv->Rho[0], dv->Rho[1], dv->Rho[2]);
    printf("Kappa=[%g %g %g] [W/(m*K)]\n", dv->Kappa[0], dv->Kappa[1], dv->Kappa[2]);
    printf("CP=[%g %g %g] [J/(kg*K)]\n", dv->CP[0], dv->CP[1], dv->CP[2]);
    printf("Mu=[%g %g %g %g] [Pa*s]\n", dv->Mu[0], dv->Mu[1], dv->Mu[2], dv->Mu[3]);

    printf("%s============ Non - dimensional params ============\n%s", KRED, KNRM);
    printf("domain_size=%g ", ndv->domain_size);
    printf("T_solid=%g\n", ndv->T_solid);
    printf("sigma_ndim=%g\n", ndv->sigma_ndim);
    printf("Grav_ndim=[%g %g %g]\n", ndv->Ggrav_ndim[0], ndv->Ggrav_ndim[1], ndv->Ggrav_ndim[2]);
    printf("Htr=%g ", ndv->Htr);
    printf("Arrhenius_const=%g ", ndv->Arrhenius_const);
    printf("Ea_by_R=%g ", ndv->Ea_by_R);
    printf("n_degree=%g\n", ndv->n_degree);
    printf("Eeta_by_Rg=%g ", ndv->Eeta_by_Rg);
    printf("chi=%g\n", ndv->chi);
    printf("rho=[%g %g %g]\n", ndv->rho[0], ndv->rho[1], ndv->rho[2]);
    printf("kappa=[%g %g %g]\n", ndv->kappa[0], ndv->kappa[1], ndv->kappa[2]);
    printf("Cp=[%g %g %g]\n", ndv->Cp[0], ndv->Cp[1], ndv->Cp[2]);
    printf("mu=[%g %g %g %g]\n", ndv->mu[0], ndv->mu[1], ndv->mu[2], ndv->mu[3]);
    printf("ratio_front_x=%g ratio_dist_x=%g ratio_dist_y=%g\n", ndv->ratio_front_x, ndv->ratio_dist_x, ndv->ratio_dist_y);
    printf("ratio_Rbmin=%g ratio_Rbmax=%g\n", ndv->ratio_Rbmin, ndv->ratio_Rbmax);
    printf("shift_x=%g shift_y=%g\n", ndv->shift_x, ndv->shift_y);
    printf("dev_r=%g develx=%g devely=%g\n", ndv->dev_r, ndv->develx, ndv->devely);
    printf("Nb=%d Ncx=%d Ncy=%d non_saturated=%d\n", ndv->Nb, ndv->Ncx, ndv->Ncy, ndv->non_saturated);

    printf("%s============ Numerical simulation params ============\n%s", KRED, KNRM);
    printf("minlevel=%d maxlevel=%d\n", numpar->minlevel, numpar->maxlevel);
    printf("iter_fp=%d\n", numpar->iter_fp);
    printf("snapshot_i=%d snapshot_t=%g dt_vtk=%g\n", numpar->snapshot_i, numpar->snapshot_t, numpar->dt_vtk);
    printf("N_smooth=%d\n", numpar->N_smooth);
    printf("TOLERANCE=%g TOLERANCEV=%g TOLERANCEP=%g\n", numpar->TOLERANCE, numpar->TOLERANCEV, numpar->TOLERANCEP);
    printf("NITERMIN=%d NITERMAX=%d\n", numpar->NITERMIN, numpar->NITERMAX);
    printf("CFL=%g CFL_ARR=%g\n", numpar->CFL, numpar->CFL_ARR);
    printf("DT=%g maxDT=%g\n", numpar->DT, numpar->maxDT);
    printf("m_bp=%g m_bp_T=%g\n", numpar->m_bp, numpar->m_bp_T);
    printf("feps=%g ueps=%g rhoeps=%g Teps=%g aeps=%g mueps=%g\n", numpar->feps, numpar->ueps, numpar->rhoeps, numpar->Teps, numpar->aeps, numpar->mueps);

    /* Free the data */
    cyaml_free(&config, &top_schema, input, 0);

    return EXIT_SUCCESS;
}