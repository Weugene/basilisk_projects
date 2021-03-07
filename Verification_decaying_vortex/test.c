// test solenoidal reconstruction of face velocity field
// i.e. refine_face_solenoidal()
scalar divutmp[];
#include "poisson.h"
#include "utils.h"
#include "../src_local/output_vtu_foreach.h"
double divmax (face vector uf)
{
  double dmax = 0.;
  foreach() {
    double d = 0.;
    foreach_dimension()
      d += uf.x[1] - uf.x[];
    // fprintf (stderr, "%g %g %g %.6f\n", x, y, z, d);
    if (fabs(d) > dmax)
      dmax = fabs(d);
  }
  return dmax;
}

int main()
{
  init_grid(64);

  origin (-0.5,-0.5,-0.5);
  foreach_dimension()
    periodic (right);

  face vector uf[];
#if dimension == 2
  // Taylor-Green vortices
  foreach_face(x)
    uf.x[] = - cos(2.*pi*x)*sin(2.*pi*y);
  foreach_face(y)
    uf.y[] =   sin(2.*pi*x)*cos(2.*pi*y);
#else // dimension == 3
  // random flow
  foreach_face()
    uf.x[] = 1. - 2.*rand()/(double)RAND_MAX;
#endif
  boundary ((scalar *){uf});
  scalar p[];
  t = 0;
  char subname[80]; sprintf(subname, "test");
  scalar l[], omega[]; foreach() l[] = level;
  output_vtu_MPI( (scalar *) {omega, p, l, divutmp}, (vector *) {uf}, subname, 0 );
#if 1
  TOLERANCE = 1e-10;
  foreach()
    p[] = 0.;
  boundary ({p});
  project (uf, p);
#endif
  t = 1;
  foreach() l[] = level;
  output_vtu_MPI( (scalar *) {omega, p, l, divutmp}, (vector *) {uf}, subname, 0 );
  fprintf (stderr, "div max before: %.9f\n", divmax(uf));
    
  uf.x.refine = refine_face_solenoidal;
  t = 2;
  foreach() l[] = level;
  output_vtu_MPI( (scalar *) {omega, p, l, divutmp}, (vector *) {uf}, subname, 0 );
  fprintf (stderr, "div max before: %.9f\n", divmax(uf));
  refine (x*x + y*y + z*z < sq(0.25) && level < 8);

#if 0  
  output_gfs (stdout);
  foreach_face(x)
    fprintf (stderr, "%g %g %g %g 0\n", x, y, z, uf.x[]);
#endif
  t = 3;
  foreach() l[] = level;
  output_vtu_MPI( (scalar *) {omega, p, l, divutmp}, (vector *) {uf}, subname, 0 );
  fprintf (stderr, "div max after: %.9f\n", divmax(uf));
}
