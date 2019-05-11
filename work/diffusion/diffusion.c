//#include "diffusion.h"
//#include "run.h"
//
//#define EPS 0.1
#define LEVEL 8
#define MINLEVEL 4
#define MAXLEVEL 10
//scalar f[];
//
//const face vector D[] = { 1. , 1. };
//
//double solution (double x, double y, double t)
//{
//    return  exp( -1. * (sq(x) + sq(y))/ (4. * t))/ (4. * 	pi* t) ;
//}
//event init (t = 0)
//{
//  	foreach()
//    	f[] = solution(x,y,0.1);
//  	boundary ({f});
//}
//event running ( i++ )
//{
//    dt = dtnext (0.01);
//    diffusion (f,dt,D);
//    boundary ({f});
//}
//event movie (t += 0.1; t <= 10) {
//    scalar l[];
//    foreach() {
//        l[] = level;
//    }
//    static FILE * fpl   = popen ("/home/weugene/basilisk/src/ppm2mpeg > grid.mpg", "w");
////    output_ppm (l, fpl, linear=true);
//        output_ppm (l, fpl, min = 0, max = MAXLEVEL);
//
//    static FILE * fpmov = popen ("/home/weugene/basilisk/src/ppm2mpeg > f.mpg", "w");
//    output_ppm (f, fpmov, linear = true);
//
//
//}
//event print ( t = 0.1 ; t += 0.1 ; t <= 1. )
//{
//    double shift = 0.1 ;
//  	// For y=0
//  	for (double x = -L0/2 ; x <= L0/2; x += L0/200.)
//    {
//        printf ("%f %f %f \n", x, interpolate (f, x, 	0.0),solution(x,0.0,t+shift));
//    }
//  	printf ("\n\n");
//
////    static FILE * fp = fopen ("f.ppm", "w");
////    output_ppm (f, fp, linear = true);
////
////    scalar l[];
////    foreach()
////    l[] = level;
////    static FILE * fpl = fopen ("grid.ppm", "w");
////    output_ppm (l, fpl, min = 0, max = MAXLEVEL);
//}
//
//
//#if TREE
//event adapt (i++) {
//  adapt_wavelet ({A}, (double[]){1e-3}, LEVEL + 1);
//}
//#endif
//int main() {
//  	// Lenght
//  	L0 = 10.;
//  	// coordinates of lower-left corner
//  	X0 = Y0 = -L0/2;
//  	//
//  	N = 128*2 ;
//
//  	run();
//}

//#include "grid/cartesian.h"
#include "run.h"
#include "diffusion.h"

scalar A[];
scalar dA[];

double dt;

int main() {
    L0 = 1.;
//    origin (-1, -1, -1);
    size (2.); //dimension of the problem!
    init_grid (1 << LEVEL);

    run();
}
A[right]  = neumann(0);
A[left]   = neumann(0);
A[top]    = neumann(0);
A[bottom] = neumann(0);

event init (t = 0) {
    foreach()
    A[] =  1./0.1*(fabs(x*x+y*y)<0.05);
    boundary ({A});
}

event integration (i++) {

    diffusion(A,dt);
    boundary ({A});

}

//event print (i=10) {
//
//    for (double y = 0 ; y <= L0; y += 0.01){
//        printf("%g %g \n",
//               y, interpolate (A, 0, y));
//    }
//
//}
event movie (t += 0.1; t <= 1) {
    scalar l[];
    foreach() {
        l[] = level;
    }
    static FILE * fpl   = popen ("/home/weugene/basilisk/src/ppm2gif > grid.gif", "w");
    //    output_ppm (l, fpl, linear=true);
    output_ppm (l, fpl, min = 0, max = MAXLEVEL);

    static FILE * fpmov = popen ("/home/weugene/basilisk/src/ppm2gif > A.gif", "w");
    output_ppm (A, fpmov, linear = false);


}
#if TREE
event adapt (i++) {
  adapt_wavelet ({A}, (double[]){1e-3}, MAXLEVEL);
}
#endif