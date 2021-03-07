#include "embed.h"
#include "navier-stokes/centered.h"

void start_screen (void);
void stop_screen (void);
char getch (void);

face vector muv[];

int main() {
  start_screen();
  L0 = 8;
  origin (-0.5, -L0/2.);
  N = 512;
  mu = muv;
  run();
  stop_screen();
}

event properties (i++) {
  foreach_face()
    muv.x[] = fm.x[]*0.125/160.;
  boundary ((scalar*){muv});
}

double U = 1;
u.n[left]  = dirichlet(U);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

u.n[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);
u.t[embed] = fabs(y) > 0.25 ? neumann(0.) : dirichlet(0.);

event init (t = 0) {
  vertex scalar phi[];
  foreach_vertex() {
    phi[] = intersection (0.5 - y, 0.5 + y);
    phi[] = intersection (phi[], sq(x) + sq(y) - sq(0.125/2.));
  }
  boundary ({phi});
  fractions (phi, cs, fs);
  foreach()
    u.x[] = cs[] ? 1 : 0;
}

event logfile (i += 5)
  printf ("%g %d %g\n", t, i, U);

event movies (t += 0.05; t <= 15.)
{
  scalar omega[], m[];
  vorticity (u, omega);
  foreach()
    m[] = cs[] - 0.5;
  boundary ({m});
  output_ppm (omega, file = "vort.mp4", box = {{-0.5,-0.5},{7.5,0.5}},
	      min = -10, max = 10, linear = true, mask = m);
}

event adapt (i++) {
  adapt_wavelet ({cs, u.x, u.y}, (double[]){1e-2, 0.03, 0.03}, 9, 4);
}

event tune_velocity (i++) {
  char a = getch();
  if (a == 'f')    //Faster
    U *= 1.05;
  if (a == 's')    //slower
    U /= 1.05;
}
