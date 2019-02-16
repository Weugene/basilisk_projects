/* interpolation on halos  */

scalar h[];

int main (int argc, char ** argv)
{
  int n = 2048;
  init_grid (n);

  double R0 = 0.1;
  foreach()
    h[] = exp(-(x*x + y*y)/(R0*R0));
  boundary ({h});
  
  /* initial coarsening (see halo.c) */
  double tolerance = 1e-4;
  while (adapt_wavelet ({h}, &tolerance, 11).nc);

#if 0
  // we reinitialise h just to make sure that trash() does its job
  trash ({h});
  foreach()
    h[] = exp(-(x*x + y*y)/(R0*R0));
  boundary ({h});
#endif

  double max = 0.;
  for (int l = 0; l < depth(); l++)
    foreach_halo (prolongation, l) 
      foreach_child() {
        double e = exp(-(x*x+y*y)/(R0*R0)) - h[];
	printf ("%g %g %d %d %g %g\n", x, y, level, cell.neighbors, h[], e);
	if (fabs(e) > max)
	  max = fabs(e);
      }

  fprintf (stderr, "maximum error on halos: %g\n", max);

  return (max > tolerance);
}
