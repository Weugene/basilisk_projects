/**
# Couette flow in a periodic channel using embeded boundaries. */
scalar divutmp[];
#define PRINT_ALL_VALUES
#include "embed.h"
#include "navier-stokes/centered.h"
//#include "view.h"
#include "../src_local/output_vtu_foreach.h"
/**
## Channel geometry */

#define width 0.5
#define EPS 1e-14

void channel (scalar cs, face vector fs) {

  vertex scalar phi[];

  /**
  The variable *EPS* is necessary to ensure that the embeded boundary
  does not coincide with a mesh interface, resulting in a division by
  a zero surface in the embeded algorithm. */
  
  foreach_vertex()
    phi[] = difference (y - L0/2. + width/2. - EPS, y - L0/2. - width/2. + EPS);
  
  boundary ({phi});
  fractions (phi, cs, fs);
}

/**
## Couette solution */

double couette (double y) {
  return ((fabs(y - L0/2.) < width/2.) ? (y-L0/2.)/(width/2.) : 0.);
}

/**
## Code */

int main() {

  /** 
  We define a left-right periodic computational fluid domain. */

//  periodic (right);
  size (1.);
  origin (0., 0.);

  /**
  We solve the Stokes equations. The value of *TOLERANCE* dictates the
  precision of the solution. */

  stokes = true;
  TOLERANCE = 1e-5;
  DT = 1.e-2;
  
//  for (N = 8; N <= 64; N *= 2) {
    N=64;
    init_grid (N);
    run ();
//  }
}

/**
The boundary condition is zero velocity on the embedded boundaries. */

u.n[embed] = dirichlet(y > L0/2. ? 1. : -1.);
u.t[embed] = dirichlet(0.);

/**
We define a reference field to establish the converence of the
solution. */

scalar un[];

/**
Initial conditions. */

event init (t = 0) {

  foreach() {
    u.x[] = 0.;
    u.y[] = 0.;
  }

  channel (cs, fs);

  /**
  The gravity vector is aligned with the channel and the viscosity is
  unity. */
  
  mu = fm;

  /**
  We use "third-order" [face flux interpolation](/src/embed.h). */

  for (scalar s in {u})
    s.third = true;
    
  foreach()
    un[] = u.x[];
}

/**
We check if the solution is stationnary. */

event steady (t += 0.1; i <= 1000) {

  double du = change (u.x, un);
 
  if (i > 1 && du < 1e-12) {
    event("vtk_file");
    return 1; /* stop */
  }
}

//event image (t = end) {
//
//  view (fov = 24, quat = {0,0,0,1}, tx = -L0/2., ty = -L0/2., bg = {1,1,1});
//  cells ();
//  draw_vof ("cs", "fs", lw = 5);
//  squares ("u.x", linear = true, spread = -1);
//  save ("ux.png");
//}

event profile (t = end) {

  static FILE * fp = fopen ("fields", "w");
  fprintf (fp, "\n\n");
  foreach()
    fprintf (fp, "%g %g %g %g %g\n", x, y, u.x[], u.y[], p[]);

  static FILE * fe = fopen ("errors", "w");
  scalar e[];
  foreach()
    e[] = u.x[] - couette (y);
  norm n = normf (e);
  fprintf (fe, "%d %g %g %g\n", N, n.avg, n.rms, n.max);
}

//Output
event vtk_file (i++){
//event vtk_file (i+=10){
//event vtk_file (t += 1){
  char subname[80]; sprintf(subname, "couette_centered");
  scalar l[], divutmpAfter[]; foreach() {l[] = level;}

  #ifdef DEBUG_MODE
  double d = 0;
      foreach() {
          d = 0.;
          foreach_dimension(){
              d += uf.x[1] - uf.x[];
          }
          divutmpAfter[] = d/(dt*Delta);
      }
      boundary((scalar*){divutmpAfter});
  #endif
  output_vtu_MPI( (scalar *) {cs, p, pf, l, divutmpAfter, divutmp},
  (vector *) {u, g},
  (vector *) {uf}, subname, t + dt );
}
/**
# Plots 

![*u.x* velocity field](example-embed-straight-channel-couette/ux.png)

~~~gnuplot Velocity profiles at all *x* positions

reset

mycolors = "dark-blue red sea-green dark-violet orange black"
mypoints = '1 2 3 4 6 7'

width = 0.5
f(x) = (x-0.5)/(width/2.)

set output 'profiles.png'
set xlabel 'u.x'
set ylabel 'y'

plot 'fields' i 0 every 5 u 3:2 ps 3 lw 2 lc rgb word(mycolors, 1) t 'N = 8', \
     'fields' i 1 every 5 u 3:2 ps 3 lw 2 lc rgb word(mycolors, 2) t 'N = 16', \
     'fields' i 2 every 5 u 3:2 ps 3 lw 2 lc rgb word(mycolors, 3) t 'N = 32', \
     'fields' i 3 every 5 u 3:2 ps 3 lw 2 lc rgb word(mycolors, 4) t 'N = 64', \
     'fields' i 0 every 5 u (f($2)):2 lc rgb 'black' t 'Couette'
~~~

~~~gnuplot Comparison of the evolution of the error norms for the flow rate *Q* with the number of cells *N*

set output 'errors.png'
set xlabel 'N'
set ylabel 'error norms'

set logscale
set xtics 8,2,64
set format y "%g"
set cbrange[1e-10:1e10]

ftitle(a,b) = sprintf('Order %4.2f', -b)
f1(x) = a1 + b1*x
f2(x) = a2 + b2*x
f3(x) = a3 + b3*x

fit f1(x) 'errors' u (log($1)):(log($2)) via a1, b1
fit f2(x) 'errors' u (log($1)):(log($3)) via a2, b2
fit f3(x) 'errors' u (log($1)):(log($4)) via a3, b3

plot "errors" u 1:2 w p lt word(mypoints, 1) lc rgb word(mycolors, 1) ps 3 lw 2 t '|Q|_1, '.ftitle(a1, b1), \
     exp (f1(log(x))) ls 1 lc rgb word(mycolors, 1) notitle, \
     "errors" u 1:3 w p lt word(mypoints, 2) lc rgb word(mycolors, 2) ps 3 lw 2 t '|Q|_2, '.ftitle(a2, b2), \
     exp (f2(log(x))) ls 1 lc rgb word(mycolors, 2) notitle, \
     "errors" u 1:4 w p lt word(mypoints, 3) lc rgb word(mycolors, 3) ps 3 lw 2 t '|Q|_{max}, '.ftitle(a3, b3), \
     exp (f3(log(x))) ls 1 lc rgb word(mycolors, 3) notitle
~~~
*/
