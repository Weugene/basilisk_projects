//#include <stdio.h>
//    #include <math.h>
//
//int main ()
//{
//  int i,j,k;
//  int N = 256 ;
//  double L = 1. ;
//  double x,y;
//  double dx = L/N ;
//  double dt = 0.00001;
//  double A[N][N];
//  double dA[N][N];
//
//
//// boundary conditions
//for (i = 0 ; i < N ; i++) A[i][0] = A[i][N-1] = 0. ;
//for (j = 0 ; j < N ; j++) A[0][j] = A[N-1][j] = 0. ;
//
//// initial conditions
//
//  for (i = 0 ; i < N ; i++)
//    {
//      for (j = 0 ; j < N ; j++)
//	{
//	  x = i*dx - 0.5 ;
//	  y = j*dx - 0.5 ;
//	  A[i][j] = 1./0.1*((fabs(x*x+y*y) < 0.05)) ;
//	}
//    }
//
// for (j = 0 ; j < N ; j++)
//	{
//	  printf("%f \n",A[(int)N/2][j]);
//	}
// printf("\n\n");
//
//  // time integration
//
//  for (k = 0 ; k < 10 ; k++)
//    {
//
//  for (i = 1 ; i < N-1 ; i++)
//    {
//       for (j = 1 ; j < N-1 ; j++)
//	{
//      dA[i][j] = (A[i+1][j] + A[i-1][j] - 2. * A[i][j])/dx/dx +
//	(A[i][j+1] + A[i][j-1] - 2. * A[i][j])/dx/dx ;
//    }
//    }
//
//  // update
//
// for (i = 0 ; i < N ; i++)
//    {
//       for (j = 0 ; j < N ; j++)
//	{
//	  A[i][j] =  A[i][j] + dt* dA[i][j] ;
//	}
//    }
//
//    }
//
//  // print solution (centerline)
//
//       for (j = 0 ; j < N ; j++)
//	{
//	  	  printf("%f \n",A[(int)N/2][j]);
//	}
//
//}


#include "grid/cartesian.h"
#include "run.h"

scalar A[];
scalar dA[];

double dt;

int main() {
    L0 = 1.;
    N = 512;


    run();
}

// boundary conditions
//bid circle;
//mask (sq(x - 0.5) + sq(y - 0.5) < sq(0.9) ? circle : none);
//A[circle] = dirichlet(0.1);
A[left] = 0.0 ;
A[top] = 0.0 ;
A[right] = 0.0 ;
A[bottom] = 0.0 ;

// initial conditions
event init (t = 0) {
    dt =0.000001;
    foreach(){
        A[]=0.0;
        A[] =  1.*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)<0.1);
    }

    boundary ({A});
}

// integration

event integration (i++) {
    foreach()
    {
        dA[] = (A[1, 0] + A[-1, 0] - 2. * A[]) / Delta / Delta +
               (A[0, 1] + A[0, -1] - 2. * A[]) / Delta / Delta;
        dt=0.1*Delta*Delta;

    }

    foreach()
    A[] = A[] + dt*dA[];
    boundary ({A});
}


//event print (t=end) {
//
//    for (double y = 0 ; y <= L0; y += 0.01){
//        printf("%g %g \n",
//               y, interpolate (A, 0, y));
//    }
////    foreach(){
////        printf("%g %g \n", y, A[]);
////    }
//
//}

event images (t+= 1.) {
//    output_ppm (A, linear = true);

    static FILE * fprho = fopen ("out.ppm", "w");
    output_ppm (A, fprho, min = 0, max = 1);
    for (double y = 0 ; y <= L0; y += 0.01){
        printf("%g %g \n",
               y, interpolate (A, 0.5, y));
    }
}
//event image (t= 0) {
//    output_ppm (A, linear = true);
//foreach() printf("A=%g",A[]);
//    static FILE * fprho = fopen ("out.ppm", "w");
//    output_ppm (A, fprho, min = 0, max = 1);
//}
//event output (t = 5) {
//    char name[80];
//    sprintf(name, "rho.dat");
//    FILE *fp = fopen(name, "w");
//    output_field({A}, fp, linear = true);
//    fclose(fp);
//}
//event movie (t += 0.2; t <= 10) {
//    static FILE * fp = popen ("ppm2mpeg > vort.mpg", "w");
//    output_ppm (A, fp, linear = true);
//}
event end (t = 10) {
    //printf ("i = %d t = %g\n", i, t);
}
