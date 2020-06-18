#include <stdio.h>
#include <math.h>

#if dimension == 2
scalar fs[];
#define ff fs[]
#else
double ff;
#endif
(const) face vector Us = {0,0,0};
(const) face vector mu = zerof,
void func(vector utau, face vector us, vector n){
    if (is_constant(Us.x)) fprintf(ferr, "!!!us is constant\n");
    if (!is_constant(Us.x)) fprintf(ferr, "!!!us is  NOT constant\n");
    foreach(){
        fprintf(ferr, "in func s=%s solid_U = %g utau=%g n=%g\n", us.x.name, us.x[], utau.x[], n.x[]);
    }
}

int main(void)
{
#if dimension > 1 && dimension < 3//sdfsdfsf
    printf("dim is %d\n", dimension)
#endif//sdfsdfsf
    ;
    L0 = 4.;
    origin (-L0/2., -L0/2.);

    N = 4;
    init_grid(N);
    face vector utau[], uf[], n[];
    foreach_face() {
        uf.x[] = 1;
        n.x[]  = x/sqrt(sq(x) + sq(y));
    }
    double ubyn = 0;
    foreach_face() {
        ubyn = 0; foreach_dimension() ubyn += uf.x[]*n.x[];
        utau.x[] = uf.x[] - ubyn;
        fprintf(ferr, "ubyn = %g utau = %g\n", ubyn, utau.x[]);
    }
    fprintf(ferr, "\n+++++++\n");
//    face vector solid_U = {1, 0};
//    foreach_face() {
//        ubyn = 0; foreach_dimension() ubyn += (uf.x[] - solid_U.x[])*n.x[];
////        utau.x[] = uf.x[] - solid_U.x[] - ubyn;
//        fprintf(ferr, "s=%s ubyn = %g utau = %g\n", utau.x.name, ubyn, utau.x[]);
//    }

    func(utau, Us, n);
    return 0;
}