#include "config_reader.h"

int main(int argc, char *argv[]){
    struct input_yaml *input;
    read_config(argc, argv, input);
}