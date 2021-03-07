/**
# Axisymmetric coordinates

For problems with a symmetry of revolution around the $z$-axis of a
[cylindrical coordinate
system](http://en.wikipedia.org/wiki/Cylindrical_coordinate_system). The
longitudinal coordinate ($z$-axis) is *x* and the radial coordinate
($\rho$- or $r$-axis) is *y*. Note that *y* (and so *Y0*) cannot be
negative.

We first define a macro which will be used in some geometry-specific
code (e.g. [curvature computation](curvature.h)). */

#define AXI 1

/**
On trees we need refinement functions. */

#if TREE
static void refine_cm_axi (Point point, scalar cm)
{
  fine(cm,0,0) = fine(cm,1,0) = y - Delta/4.;
  fine(cm,0,1) = fine(cm,1,1) = y + Delta/4.;
}

static void refine_face_x_axi (Point point, scalar fm)
{
  if (!is_refined(neighbor(-1))) {
    fine(fm,0,0) = y - Delta/4.;
    fine(fm,0,1) = y + Delta/4.;
  }
  if (!is_refined(neighbor(1)) && neighbor(1).neighbors) {
    fine(fm,2,0) = y - Delta/4.;
    fine(fm,2,1) = y + Delta/4.;
  }
  fine(fm,1,0) = y - Delta/4.;
  fine(fm,1,1) = y + Delta/4.;
}

static void refine_face_y_axi (Point point, scalar fm)
{
  if (!is_refined(neighbor(0,-1)))
    fine(fm,0,0) = fine(fm,1,0) = max(y - Delta/2., 1e-20);
  if (!is_refined(neighbor(0,1)) && neighbor(0,1).neighbors)
    fine(fm,0,2) = fine(fm,1,2) = y + Delta/2.;
  fine(fm,0,1) = fine(fm,1,1) = y;
}
#endif

event metric (i = 0) {

  /**
  By default *cm* is a constant scalar field. To make it variable, we
  need to allocate a new field. We also move it at the begining of the
  list of variables: this is important to ensure the metric is defined
  before other fields. */

  if (is_constant(cm)) {
    scalar * l = list_copy (all);
    cm = new scalar;
    free (all);
    all = list_concat ({cm}, l);
    free (l);
  }

  /**
  The volume/area of a cell is proportional to $r$ (i.e. $y$). We need
  to set boundary conditions at the top and bottom so that *cm* is
  interpolated properly when refining/coarsening the mesh. */

  scalar cmv = cm;
  foreach()
    cmv[] = y;
  cm[top] = dirichlet(y);
  cm[bottom] = dirichlet(y);

  /**
  We do the same for the length scale factors. The "length" of faces
  on the axis of revolution is zero ($y=r=0$ on the axis). To avoid
  division by zero we set it to epsilon (note that mathematically the
  limit is well posed). */

  if (is_constant(fm.x)) {
    scalar * l = list_copy (all);
    fm = new face vector;
    free (all);
    all = list_concat ((scalar *){fm}, l);
    free (l);
  }
  face vector fmv = fm;
  foreach_face()
    fmv.x[] = max(y, 1e-20);
  fm.t[top] = dirichlet(y);
  fm.t[bottom] = dirichlet(y);
  
  /**
  We set our refinement/prolongation functions on trees. */

#if TREE
  cm.refine = cm.prolongation = refine_cm_axi;
  fm.x.prolongation = refine_face_x_axi;
  fm.y.prolongation = refine_face_y_axi;
#endif
  
  boundary ({cm, fm});
}

/**
## See also

* [Axisymmetric streamfunction](axistream.h)
*/
