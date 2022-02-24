/******************************************************************************
* CYAML schema to tell libcyaml about both expected YAML and data structure.
*
* (Our CYAML schema is just a bunch of static const data.)
******************************************************************************/

//#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <cyaml/cyaml.h>
// nondimensional numbers
struct numbers {
    double Re; // Reynolds number Re = rho*U*d/mu
    double Ca; // Capillary number Ca = Mu*Ud/sigma
};

// dimensional variables
struct dimensional_vars {
    double characteristic_size;
};

struct input_yaml {
    char* name;
    struct numbers n;
    struct dimensional_vars dv;
};

//struct test = {
//    float foo;
//    int bar[10];
//}

//static const cyaml_schema_value_t bar_schema = {
//        CYAML_VALUE_INT(CYAML_FLAG_DEFAULT, int)
//};

//static const cyaml_schema_field_t test_fields_schema[] = {
//        CYAML_FIELD_FLOAT("name", CYAML_FLAG_OPTIONAL, struct test, foo),
//        CYAML_FIELD_FLOAT("n", CYAML_FLAG_OPTIONAL, struct numbers, n),
////        CYAML_FIELD_SEQUENCE_FIXED("bar", CYAML_FLAG_DEFAULT, struct test, bar, &bar_schema, 10),
//        CYAML_FIELD_END
//};

//static const cyaml_schema_value_t test_schema = {
//        CYAML_VALUE_MAPPING(CYAML_FLAG_POINTER, struct test, test_fields_schema)
//};
/* CYAML value schema for entries of the data sequence. */
//static const cyaml_schema_value_t data_entry = {
//        CYAML_VALUE_FLOAT(CYAML_FLAG_DEFAULT, double),
//};


static const cyaml_schema_field_t numbers_fields[] = {
        CYAML_FIELD_FLOAT("Re", CYAML_FLAG_DEFAULT, struct numbers, Re),
        CYAML_FIELD_FLOAT("Ca", CYAML_FLAG_DEFAULT, struct numbers, Ca),
};

static const cyaml_schema_field_t dimensional_fields[] = {
        CYAML_FIELD_FLOAT("characteristic_size", CYAML_FLAG_DEFAULT, struct dimensional_vars, characteristic_size),
};

/* CYAML mapping schema fields array for the top level mapping. */
static const cyaml_schema_field_t top_mapping_schema[] = {
        CYAML_FIELD_STRING_PTR("name", CYAML_FLAG_POINTER,struct input_yaml, name,0, CYAML_UNLIMITED),
        CYAML_FIELD_MAPPING("n", CYAML_FLAG_DEFAULT,struct input_yaml, n, numbers_fields),
//        CYAML_FIELD_MAPPING("dv", CYAML_FLAG_DEFAULT,struct input_yaml, dv, dimensional_fields),
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

int main(int argc, char *argv[])
{
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
//    printf("%s:\n", n->name);
//    for (unsigned i = 0; i < n->data_count; i++) {
//        printf("  - %i\n", n->data[i]);
//    }

    /* Free the data */
    cyaml_free(&config, &top_schema, n, 0);

    return EXIT_SUCCESS;
}