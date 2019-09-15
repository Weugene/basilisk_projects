#include "run.h"


int main() {
    L0 = 8.;
    origin (-L0/2, -L0/2.);
    N = 64;
    run();
}

event init (t = 0) {
printf("0**1.6=%g\n", pow(1e-30,1.6));
}

event vtk_file (i=6; i<=10; i+=1)
{
    printf("Here output only after i=%d>=6\n", i)
}

event vtk_file_other (i+=1; i<6)
{
    printf("Here output i=%d<6\n", i)
}


event stop(i = 10);
