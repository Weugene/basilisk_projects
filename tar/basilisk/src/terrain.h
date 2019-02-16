#include <stdarg.h>
#include <kdt/kdt.h>
#pragma autolink -L$BASILISK/kdt -lkdt

@if _OPENMP
@ define NPROC omp_get_max_threads()
@ define PROC tid()
@else
@ define NPROC 1
@ define PROC 0
@endif

attribute {
  void ** kdt;
  scalar n, dmin, dmax;
}

static int includes (KdtRect rect, Point * p)
{
  Point point = *p;
  Delta_x /= 2.; Delta_y /= 2.;
  return (rect[0].l >= x - Delta_x && rect[0].h <= x + Delta_x &&
	  rect[1].l >= y - Delta_y && rect[1].h <= y + Delta_y);
}

static int intersects (KdtRect rect, Point * p)
{
  Point point = *p;
  Delta_x /= 2.; Delta_y /= 2.;
  return (rect[0].l <= x + Delta_x && rect[0].h >= x - Delta_x &&
	  rect[1].l <= y + Delta_y && rect[1].h >= y - Delta_y);
}

static void reconstruct_terrain (Point point, scalar zb)
{
  KdtSum s;
  kdt_sum_init (&s);
  Delta_x /= 2.; Delta_y /= 2.;
  KdtRect rect = {{x - Delta_x, x + Delta_x},
		  {y - Delta_y, y + Delta_y}};
  for (Kdt ** kdt = (Kdt **) zb.kdt[PROC]; *kdt; kdt++)
    kdt_query_sum (*kdt,
		   (KdtCheck) includes,
		   (KdtCheck) intersects, &point,
		   rect, &s);
  scalar n = zb.n, dmin = zb.dmin, dmax = zb.dmax;
  n[] = s.n;
  if (s.w > 0.) {
    zb[] = s.H0/s.w;
    dmin[] = s.Hmin;
    dmax[] = s.Hmax;
  }
  else {
    /* not enough points in database, use bilinear interpolation
       from coarser level instead */
    if (level > 0)
      zb[] = (9.*coarse(zb,0,0) + 
	      3.*(coarse(zb,child.x,0) + coarse(zb,0,child.y)) + 
	      coarse(zb,child.x,child.y))/16.;
    else
      zb[] = 0.; // no points at level 0!
    dmin[] = nodata;
    dmax[] = nodata;
  }
}

void refine_terrain (Point point, scalar zb)
{
  foreach_child()
    reconstruct_terrain (point, zb);
}

static void delete_terrain (scalar zb)
{
  for (int i = 0; i < NPROC; i++) {
    for (Kdt ** kdt = (Kdt **) zb.kdt[i]; *kdt; kdt++)
      kdt_destroy (*kdt);
    free (zb.kdt[i]);
  }
  free (zb.kdt);
}

void terrain (scalar zb, ...)
{
  zb.kdt = qcalloc (NPROC, void *);
  zb.delete = delete_terrain;
  
  int nt = 0;
  va_list ap;
  va_start (ap, zb);
  char * name;
  while ((name = va_arg (ap, char *))) {
    for (int i = 0; i < NPROC; i++) {
      Kdt ** kdt = (Kdt **) zb.kdt[i];
      zb.kdt[i] = qrealloc (kdt, nt + 2, Kdt *);
      kdt[nt] = kdt_new();
      kdt[nt + 1] = NULL;
      if (kdt_open (kdt[nt], name)) {
	fprintf (stderr, "terrain: could not open terrain database '%s'\n", 
		 name);
	exit (1);
      }
    }
    nt++;
  }
  va_end (ap);

  scalar n = new scalar;
  scalar dmin = new scalar;
  scalar dmax = new scalar;
  zb.n = n;
  zb.dmin = dmin;
  zb.dmax = dmax;

#if TREE
  zb.refine = refine_terrain;
  n.refine = no_restriction;
  n.restriction = no_restriction;
  dmin.refine = no_restriction;
  dmin.restriction = no_restriction;
  dmax.refine = no_restriction;
  dmax.restriction = no_restriction;
#endif

  trash ({zb});
  for (int l = 0; l <= depth(); l++) {
    foreach_level (l)
      reconstruct_terrain (point, zb);
    boundary_level ({zb}, l);
  }
}
