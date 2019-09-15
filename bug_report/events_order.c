#include "run.h"


int main() {
    L0 = 8.;
    origin (-L0/2, -L0/2.);
    N = 64;
    run();
}

event init (i = 0) {
    fprintf(stderr, "0**1.6=%g\n", pow(1e-30,1.6));
}

event example (i+=1)
{
    fprintf(stderr, "example 1! i=%d\n", i);
}

event example (i+=1,last)
{
    fprintf(stderr, "example 2! i=%d\n", i);
}

event example ( i+=1,last)
{
    fprintf(stderr, "example 3! i=%d\n", i);
}

event example ( i+=1)
{
    fprintf(stderr, "example 4! i=%d\n", i);
}


event exam (i+=1; i<6)
{
    fprintf(stderr, "exam! i=%d\n", i);
}

event ex (i+=1; i<6)
{
    fprintf(stderr, "ex! i=%d\n", i);
}
event stop(i = 10);

