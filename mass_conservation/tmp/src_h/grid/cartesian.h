#define GRIDNAME "Cartesian"
#define dimension 2
#define GHOSTS 1

#define _I     (point.i - 1)
#define _J     (point.j - 1)
#define _DELTA (1./point.n)

typedef struct {
  Grid g;
  char * d;
  int n;
} Cartesian;

struct _Point {
  int i, j, level, n;
};
static Point last_point;

#define cartesian ((Cartesian *)grid)

@def data(k,l,m) ((double *)&cartesian->d[((point.i + k)*(point.n + 2) +
					 (point.j + l))*datasize]) @
@define allocated(...) true

@define POINT_VARIABLES VARIABLES

@def foreach()
  OMP_PARALLEL() {
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  Point point;
  point.n = cartesian->n;
  int _k;
  OMP(omp for schedule(static))
  for (_k = 1; _k <= point.n; _k++) {
    point.i = _k;
    for (point.j = 1; point.j <= point.n; point.j++) {
      POINT_VARIABLES
@
@define end_foreach() }}}

@def foreach_face_generic()
  OMP_PARALLEL() {
  int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);
  Point point;
  point.n = cartesian->n;
  int _k;
  OMP(omp for schedule(static))
  for (_k = 1; _k <= point.n + 1; _k++) {
    point.i = _k;
    for (point.j = 1; point.j <= point.n + 1; point.j++) {
      POINT_VARIABLES
@
@define end_foreach_face_generic() }}}

@def foreach_vertex()
foreach_face_generic() {
  x -= Delta/2.; y -= Delta/2.;
@
@define end_foreach_vertex() } end_foreach_face_generic()

#define foreach_edge() foreach_face(y,x)

@define is_face_x() (point.j <= point.n)
@define is_face_y() (point.i <= point.n)

@if TRASH
@ undef trash
@ define trash(list) reset(list, undefined)
@endif

#include "neighbors.h"

void reset (void * alist, double val)
{
  scalar * list = (scalar *) alist;
  for (int i = 0; i < sq(cartesian->n + 2); i++)
    for (scalar s in list)
      if (!is_constant(s))
	((double *)(&cartesian->d[i*datasize]))[s.i] = val;
}

// Boundaries

@def foreach_boundary_dir(l,d)
  OMP_PARALLEL() {
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  Point point;
  point.n = cartesian->n;
  int * _i = &point.j;
  if (d == left) {
    point.i = GHOSTS;
    ig = -1;
  }
  else if (d == right) {
    point.i = point.n + GHOSTS - 1;
    ig = 1;
  }
  else if (d == bottom) {
    point.j = GHOSTS;
    _i = &point.i;
    jg = -1;
  }
  else if (d == top) {
    point.j = point.n + GHOSTS - 1;
    _i = &point.i;
    jg = 1;
  }
  int _l;
  OMP(omp for schedule(static))
  for (_l = 0; _l < point.n + 2*GHOSTS; _l++) {
    *_i = _l;
    {
      POINT_VARIABLES
@
@def end_foreach_boundary_dir()
    }
  }
  }
@

@define neighbor(o,p,q) ((Point){point.i+o, point.j+p, point.level, point.n})
@def is_boundary(point) (point.i < GHOSTS || point.i >= point.n + GHOSTS ||
			 point.j < GHOSTS || point.j >= point.n + GHOSTS)
@

@def foreach_boundary(b)
  foreach_boundary_dir (depth(), b)
    if (!is_boundary(point)) {
@
@define end_foreach_boundary() } end_foreach_boundary_dir()

// ghost cell coordinates for each direction
static int _ig[] = {1,-1,0,0}, _jg[] = {0,0,1,-1};

static void box_boundary_level_normal (const Boundary * b, scalar * list, int l)
{
  int d = ((BoxBoundary *)b)->d;

  OMP_PARALLEL() {
    Point point;
    point.n = cartesian->n;
    if (d % 2)
      ig = jg = 0;
    else {
      ig = _ig[d]; jg = _jg[d];
    }
    int _start = GHOSTS, _end = point.n + GHOSTS, _k;  
    OMP(omp for schedule(static))
      for (_k = _start; _k < _end; _k++) {
	point.i = d > left ? _k : d == right ? point.n + GHOSTS - 1 : GHOSTS;
	point.j = d < top  ? _k : d == top   ? point.n + GHOSTS - 1 : GHOSTS;
	Point neighbor = {point.i + ig, point.j + jg};
	for (scalar s in list) {
	  scalar b = s.v.x;
	  val(s,ig,jg) = b.boundary[d] (point, neighbor, s, NULL);
	}
      }
  }
}

static void box_boundary_level_tangent (const Boundary * b, 
					scalar * list, int l)
{
  int d = ((BoxBoundary *)b)->d;

  OMP_PARALLEL() {
    Point point;
    point.n = cartesian->n;
    ig = _ig[d]; jg = _jg[d];
    int _start = GHOSTS, _end = point.n + 2*GHOSTS, _k;
  
    OMP(omp for schedule(static))
      for (_k = _start; _k < _end; _k++) {
	point.i = d > left ? _k : d == right ? point.n + GHOSTS - 1 : GHOSTS;
	point.j = d < top  ? _k : d == top   ? point.n + GHOSTS - 1 : GHOSTS;
	Point neighbor = {point.i + ig, point.j + jg};
	for (scalar s in list) {
	  scalar b = s.v.y;
	  val(s,ig,jg) = b.boundary[d] (point, neighbor, s, NULL);
	}
      }
  }
}

static void box_boundary_level (const Boundary * b, scalar * list, int l)
{
  int d = ((BoxBoundary *)b)->d;
  scalar * centered = NULL, * normal = NULL, * tangent = NULL;

  int component = d/2;
  for (scalar s in list)
    if (!is_constant(s)) {
      if (s.face) {
	if ((&s.d.x)[component]) {
	  scalar b = s.v.x;
	  if (b.boundary[d])
	    normal = list_add (normal, s);
	}
	else {
	  scalar b = s.v.y;
	  if (b.boundary[d])
	    tangent = list_add (tangent, s);
	}
      }	
      else if (s.boundary[d])
	centered = list_add (centered, s);
    }

  OMP_PARALLEL() {
    Point point;
    point.n = cartesian->n;
    ig = _ig[d]; jg = _jg[d];
    int _start = 1, _end = point.n, _k;
    /* traverse corners only for top and bottom */
    if (d > left) { _start--; _end++; }
    OMP(omp for schedule(static))
      for (_k = _start; _k <= _end; _k++) {
	point.i = d > left ? _k : d == right ? point.n : 1;
	point.j = d < top  ? _k : d == top   ? point.n : 1;
	Point neighbor = {point.i + ig, point.j + jg};
	for (scalar s in centered) {
	  scalar b = (s.v.x.i < 0 ? s :
		      s.i == s.v.x.i && d < top ? s.v.x :
		      s.i == s.v.y.i && d >= top ? s.v.x :
		      s.v.y);
	  val(s,ig,jg) = b.boundary[d] (point, neighbor, s, NULL);
	}
      }
  }
  free (centered);

  box_boundary_level_normal (b, normal, l);
  free (normal);
  box_boundary_level_tangent (b, tangent, l);
  free (tangent);
}

void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
  free (cartesian->d);
  free (cartesian);
  grid = NULL;
}

void init_grid (int n)
{
  if (cartesian && n == cartesian->n)
    return;
  free_grid();
  Cartesian * p = qmalloc (1, Cartesian);
  size_t len = (n + 2)*(n + 2)*datasize;
  p->n = N = n;
  p->d = qmalloc (len, char);
  /* trash the data just to make sure it's either explicitly
     initialised or never touched */
  double * v = (double *) p->d;
  for (int i = 0; i < len/sizeof(double); i++)
    v[i] = undefined;
  grid = (Grid *) p;
  reset (all, 0.);
  for (int d = 0; d < nboundary; d++) {
    BoxBoundary * box = qcalloc (1, BoxBoundary);
    box->d = d;
    Boundary * b = (Boundary *) box;
    b->level   = box_boundary_level;
    add_boundary (b);
  }
  // mesh size
  grid->n = grid->tn = sq(n);
}

void realloc_scalar (int size)
{
  Cartesian * p = cartesian;
  size_t len = (p->n + 2)*(p->n + 2);
  qrealloc (p->d, len*(datasize + size), char);
  char * data = p->d + (len - 1)*datasize;
  for (int i = len - 1; i > 0; i--, data -= datasize)
    memmove (data + i*size, data, datasize);
  datasize += size;
}

struct _locate { double x, y, z; };

Point locate (struct _locate p)
{
  Point point;
  point.n = cartesian->n;
  point.i = (p.x - X0)/L0*point.n + 1;
  point.j = (p.y - Y0)/L0*point.n + 1;
  point.level = (point.i >= 1 && point.i <= point.n &&
		 point.j >= 1 && point.j <= point.n) ? 0 : - 1;
  return point;
}

#include "cartesian-common.h"
