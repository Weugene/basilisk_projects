/**
# Vertical remapping

This implements a simple vertical remapping to "$\sigma$-coordinates"
(equally-distributed by default).

We use the [PPR Library](https://github.com/dengwirda/PPR) of Engwirda
and Kelley to perform the remapping. The default settings are using
the Parabolic Piecewise Method without limiting. */

#include "ppr/ppr.h"

// int edge_meth = p1e_method, cell_meth = plm_method, cell_lim = null_limit;
int edge_meth = p3e_method, cell_meth = ppm_method, cell_lim = null_limit;
// int edge_meth = p5e_method, cell_meth = pqm_method, cell_lim = null_limit;

/**
The distribution of layers can be controlled using the *beta* array
which defines the ratio of the thickness of each layer to the total
depth $H$ (i.e. the relative thickness). By default all layers have
the same relative thickness. */

double * beta = NULL;

event defaults (i = 0)
{
  beta = malloc (nl*sizeof(double));
  for (int l = 0; l < nl; l++)
    beta[l] = 1./nl;
}

/**
The *vertical_remapping()* function takes a (block) field of layer
thicknesses and the corresponding list of tracer fields and performs
the remapping (defined by *beta*). */

trace
void vertical_remapping (scalar h, scalar * tracers)
{
  int nvar = list_len(tracers), ndof = 1, npos = nl + 1;
  foreach() {
    double H = 0.;
    foreach_layer()
      H += h[];

    if (H > dry) {
      double zpos[npos], znew[npos];
      double fdat[nvar*nl], fnew[nvar*nl];
      zpos[0] = znew[0] = 0.;
      foreach_layer() {
	zpos[point.l+1] = zpos[point.l] + h[];
	int i = nvar*point.l;
	for (scalar s in tracers)
	  fdat[i++] = s[];
	h[] = H*beta[point.l];
	znew[point.l+1] = znew[point.l] + h[];
      }

      my_remap (&npos, &npos, &nvar, &ndof, zpos, znew, fdat, fnew,
		&edge_meth, &cell_meth, &cell_lim);

      foreach_layer() {
	int i = nvar*point.l;
	for (scalar s in tracers)
	  s[] = fnew[i++];
      }
    }
  }
  boundary ({h});
  boundary (tracers);
}

/**
The remapping is applied at every timestep just before the vertical
viscosity i.e. just after horizontal advection in the [multilayer
solver](hydro.h). */

event remap (i++) {
  if (nl > 1)
    vertical_remapping (h, tracers);
}

/**
The *beta* array is freed at the end of the run. */

event cleanup (i = end)
{
  free (beta), beta = NULL;
}
