/******************************************************************************
* CYAML schema to tell libcyaml about both expected YAML and data structure.
*
* (Our CYAML schema is just a bunch of static const data.)
******************************************************************************/

#pragma once
#include <cyaml/cyaml.h>
// nondimensional numbers
struct numbers {
    double Re; // Reynolds number Re = rho*U*d/mu
    double Ca; // Capillary number Ca = Mu*Ud/sigma
    double Pe; // Peclet number Pe = rho*Cp*U*d/kappa
    double Pr; // Prandtl number
    double Fr; // Froude number Fr = sqrt(u^2/(g*cyl_diam))
};

// dimensional variables
struct dimensional_vars {
    double characteristic_size;
    double Uin;
    double Ggrav;
    double Tin;
    double T_solid;
    double Tam;
    double Rho1, Rho2, Rho3;
    double Mu0, Mu1, Mu2, Mu3;
    double Kappa1, Kappa2, Kappa3;
    double CP1, CP2, CP3;
    double Sigma;
    double Htr;
    double Arrhenius_const;
    double Ea_by_R;
    double n_degree;
    double Eeta_by_Rg;
    double chi;
};

// dimensionless variables
struct dimensionless_vars {
    double characteristic_size;
    double Uin;
    double rho1, rho2, rho3;
    double mu0, mu1, mu2, mu3;
    double kappa1, kappa2, kappa3;
    double Cp1, Cp2, Cp3;
    double sigma;
    double RhoR, RhoRS, MuR, MuRS, CpR, CpRS, KappaR, KappaRS

};

struct sim_params {
    double TOLERANCE;
    double TOLERANCEV;
    double TOLERANCEP;
    int NITERMIN;
    int NITERMAX;
    int N_smooth;
    double CFL;
    double CFL_ARR;
    double DT;
    double maxDT0;
    bool stokes;
    bool stokes_heat;
    double m_bp;
    double m_bp_T;
};

struct input_yaml {
    char* name;
    numbers n;
    dimensional_vars dv;
    dimensionless_vars ndv;
    sim_params simp;
};

/* CYAML value schema for entries of the data sequence. */
static const cyaml_schema_value_t data_entry = {
        CYAML_VALUE_FLOAT(CYAML_FLAG_DEFAULT, int),
};

static const cyaml_schema_field_t numbers_fields[] = {
        CYAML_FIELD_FLOAT(
                "Re", CYAML_FLAG_OPTIONAL,
        double, Re),
}

/* CYAML mapping schema fields array for the top level mapping. */
static const cyaml_schema_field_t top_mapping_schema[] = {
        CYAML_FIELD_STRING_PTR("name", CYAML_FLAG_POINTER,
        struct input_yaml, name,
        0, CYAML_UNLIMITED),
        CYAML_FIELD_FLOAT_PTR("Re", CYAML_FLAG_POINTER,
        struct input_yaml, name),
        CYAML_FIELD_SEQUENCE("data", CYAML_FLAG_POINTER,
        struct input_yaml, data, &data_entry,
        0, CYAML_UNLIMITED),
        CYAML_FIELD_END
};

/* CYAML value schema for the top level mapping. */
static const cyaml_schema_value_t top_schema = {
        CYAML_VALUE_MAPPING(CYAML_FLAG_POINTER,
        struct input_yaml, top_mapping_schema),
};


/******************************************************************************
 * Actual code to load and save YAML doc using libcyaml.
 ******************************************************************************/

/* Our CYAML config.
 *
 * If you want to change it between calls, don't make it const.
 *
 * Here we have a very basic config.
 */
static const cyaml_config_t config = {
        .log_fn = cyaml_log,            /* Use the default logging function. */
        .mem_fn = cyaml_mem,            /* Use the default memory allocator. */
        .log_level = CYAML_LOG_WARNING, /* Logging errors and warnings only. */
};

int read_config(int argc, char *argv[])){

    cyaml_err_t err;
    struct input_yaml *n;
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
                          &top_schema, (cyaml_data_t **)&n, NULL);
    if (err != CYAML_OK) {
        fprintf(stderr, "ERROR: %s\n", cyaml_strerror(err));
        return EXIT_FAILURE;
    }

    /* Use the data. */
    printf("%s:\n", n->name);
    for (unsigned i = 0; i < n->data_count; i++) {
        printf("  - %i\n", n->data[i]);
    }

    /* Free the data */
    cyaml_free(&config, &top_schema, n, 0);

    return EXIT_SUCCESS;
}