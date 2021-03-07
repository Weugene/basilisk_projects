#include <stdio.h>
#include <math.h>

#if dimension == 2
scalar fs[];
#define ff fs[]
#else
double ff;
#endif
int main(void)
{
    #if dimension > 1 && dimension < 3//sdfsdfsf
        printf("dim is %d\n", dimension)
    #endif;//sdfsdfsf
    ;
        L0 = 8.;
        origin (-0.5, -L0/2.);

        N = 4;
        init_grid(N);

        foreach() {
            ff = x;
            printf("dim is %g\n", ff);
        }
    return 0;
}


