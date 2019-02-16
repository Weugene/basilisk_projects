#line 1 "atomisation-cpp.c"
#line 1 "<built-in>"
#line 1 "<command-line>"
#line 1 "/usr/include/stdc-predef.h"
#line 1 "<command-line>"
#line 1 "atomisation-cpp.c"
#if _XOPEN_SOURCE < 700
#undef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif
#if _GNU_SOURCE
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#endif



#line 1 "/home/popinet/basilisk/src/common.h"
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#define pi 3.14159265358979
#undef HUGE
#define HUGE ((double)1e30)
#define nodata HUGE
#define _NVARMAX 65536
#define is_constant(v) ((v).i >= _NVARMAX)
#define constant(v) (is_constant(v) ? _constant[(v).i - _NVARMAX] : nodata)

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define sq(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))
#define sign(x) ((x) > 0 ? 1 : -1)
#define noise() (1. - 2.*rand()/(double)RAND_MAX)
#define clamp(x,a,b) ((x) < (a) ? (a) : (x) > (b) ? (b) : (x))
#define swap(type,a,b) { type tmp = a; a = b; b = tmp; }
#define unmap(x,y)

#define trash(x)


#define systderr stderr
#define systdout stdout

#if _MPI
static FILE * qstderr (void);
static FILE * qstdout (void);
FILE * ferr, * fout;
#else
# define qstderr() stderr
# define qstdout() stdout
# define ferr stderr
# define fout stdout
#endif



#if MTRACE

struct {
  FILE * fp;
  size_t total, max;
  size_t overhead, maxoverhead;
  size_t nr;
  size_t startrss, maxrss;
  char * fname;
} pmtrace;

typedef struct {
  char * func, * file;
  size_t max, total;
  int line, id;
} pmfunc;

typedef struct {
  size_t id, size;
} pmdata;

static pmfunc * pmfuncs = NULL;
static int pmfuncn = 0;

#define sysmalloc malloc
#define syscalloc calloc
#define sysrealloc realloc
#define sysfree free
#define systrdup strdup

static int pmfunc_index (const char * func, const char * file, int line)
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++)
    if (p->line == line && !strcmp(func, p->func) && !strcmp(file, p->file))
      return p->id;
  pmfuncn++;
  pmfuncs = sysrealloc (pmfuncs, pmfuncn*sizeof(pmfunc));
  p = &pmfuncs[pmfuncn - 1];
  memset (p, 0, sizeof(pmfunc));
  p->func = systrdup(func);
  p->file = systrdup(file);
  p->line = line;
  p->id = pmfuncn;
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "@ %d %s %s %d\n", pmfuncn, func, file, line);
  return pmfuncn;
}

static void pmfunc_trace (pmfunc * f, char c)
{
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "%c %d %ld %ld %ld",
      c, f->id, pmtrace.nr, pmtrace.total, f->total);
#if (_GNU_SOURCE || _DARWIN_C_SOURCE)
  if (pmtrace.nr % 1 == 0) {
    struct rusage usage;
    getrusage (RUSAGE_SELF, &usage);
    if (pmtrace.fp)
      fprintf (pmtrace.fp, " %ld", usage.ru_maxrss*1024);
    if (!pmtrace.nr)
      pmtrace.startrss = usage.ru_maxrss;
    if (usage.ru_maxrss > pmtrace.maxrss)
      pmtrace.maxrss = usage.ru_maxrss;
  }
#endif
  if (pmtrace.fp)
    fputc ('\n', pmtrace.fp);
  pmtrace.nr++;
}

static void * pmfunc_alloc (pmdata * d, size_t size,
       const char * func, const char * file, int line,
       char c)
{
  assert (d != NULL);
  d->id = pmfunc_index(func, file, line);
  d->size = size;
  pmfunc * f = &pmfuncs[d->id - 1];
  f->total += size;
  if (f->total > f->max)
    f->max = f->total;
  pmtrace.total += size;
  pmtrace.overhead += sizeof(pmdata);
  if (pmtrace.total > pmtrace.max) {
    pmtrace.max = pmtrace.total;
    pmtrace.maxoverhead = pmtrace.overhead;
  }
  pmfunc_trace (f, c);
  return ((char *)d) + sizeof(pmdata);
}

static void * pmfunc_free (void * ptr, char c)
{
  if (!ptr)
    return ptr;
  pmdata * d = (pmdata *) (((char *)ptr) - sizeof(pmdata));
  if (d->id < 1 || d->id > pmfuncn) {
    fputs ("*** MTRACE: ERROR!: corrupted free()", qstderr());
    if (d->size == 0)
      fputs (", possible double free()", qstderr());
    else
      fputs (", not traced?", qstderr());
    fputs (", aborting...\n", qstderr());
    abort();
    return ptr;
  }
  else {
    pmfunc * f = &pmfuncs[d->id - 1];
    if (f->total < d->size) {
      fprintf (qstderr(), "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        f->total, d->size);
      abort();
    }
    else
      f->total -= d->size;
    if (pmtrace.total < d->size) {
      fprintf (qstderr(), "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        pmtrace.total, d->size);
      abort();
    }
    else {
      pmtrace.total -= d->size;
      pmtrace.overhead -= sizeof(pmdata);
    }
    d->id = 0;
    d->size = 0;
    pmfunc_trace (f, c);
    return d;
  }
}

static void * pmalloc (size_t size,
         const char * func, const char * file, int line)
{
  return pmfunc_alloc (sysmalloc (sizeof(pmdata) + size),
         size, func, file, line, '+');
}

static void * pcalloc (size_t nmemb, size_t size,
         const char * func, const char * file, int line)
{
  void * p = pmalloc (nmemb*size, func, file, line);
  return memset (p, 0, nmemb*size);
}

static void * prealloc (void * ptr, size_t size,
   const char * func, const char * file, int line)
{
  return pmfunc_alloc (sysrealloc (pmfunc_free(ptr, '<'),
       sizeof(pmdata) + size),
         size, func, file, line, '>');
}

static void pfree (void * ptr,
     const char * func, const char * file, int line)
{
  sysfree (pmfunc_free (ptr, '-'));
}

static char * pstrdup (const char * s,
         const char * func, const char * file, int line)
{
  char * d = pmalloc (strlen(s) + 1, func, file, line);
  return strcpy (d, s);
}

#if MTRACE < 3
static int pmaxsort (const void * a, const void * b) {
  const pmfunc * p1 = a, * p2 = b;
  return p1->max < p2->max;
}
#endif

static int ptotalsort (const void * a, const void * b) {
  const pmfunc * p1 = a, * p2 = b;
  return p1->total < p2->total;
}

static void pmfuncs_free()
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++) {
    sysfree (p->func);
    sysfree (p->file);
  }
  sysfree (pmfuncs);
}

void pmuntrace (void)
{
#if MTRACE < 3
  fprintf (qstderr(),
    "*** MTRACE: max resident  set size: %10ld bytes\n"
    "*** MTRACE: max traced memory size: %10ld bytes"
    " (tracing overhead %.1g%%)\n"
    "%10s    %20s   %s\n",
    pmtrace.maxrss*1024,
    pmtrace.max, pmtrace.maxoverhead*100./pmtrace.max,
    "max bytes", "function", "file");
  qsort (pmfuncs, pmfuncn, sizeof(pmfunc), pmaxsort);
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn && p->max > 0; i++, p++)
    fprintf (qstderr(), "%10ld    %20s   %s:%d\n",
      p->max, p->func, p->file, p->line);

  if (pmtrace.fp) {
    char * fname = pmtrace.fname, * s;
    while ((s = strchr(fname,'/')))
      fname = s + 1;

    fputs ("load(\"`echo $BASILISK`/mtrace.plot\")\n", pmtrace.fp);
    fprintf (pmtrace.fp,
      "plot '%s' u 3:($6-%g) w l t 'ru_maxrss - %.3g',"
      "total(\"%s\") w l t 'total'",
      fname,
      pmtrace.startrss*1024.,
      pmtrace.startrss*1024.,
      fname);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->max > 0.01*pmtrace.max; i++, p++)
      fprintf (pmtrace.fp,
        ",func(\"%s\",%d) w l t '%s'",
        fname, p->id, p->func);
    fputc ('\n', pmtrace.fp);
    fprintf (qstderr(),
      "*** MTRACE: To get a graph use: tail -n 2 %s | gnuplot -persist\n",
      fname);
    fclose (pmtrace.fp);
    pmtrace.fp = NULL;
    sysfree (pmtrace.fname);
  }
#endif

  if (pmtrace.total > 0) {
    qsort (pmfuncs, pmfuncn, sizeof(pmfunc), ptotalsort);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->total > 0; i++, p++)
      fprintf (qstderr(), "%s:%d: error: %ld bytes leaked here\n",
        p->file, p->line, p->total);
    pmfuncs_free();
    exit(1);
  }
  else {
#if MTRACE < 3
    fputs ("*** MTRACE: No memory leaks\n", qstderr());
#endif
    pmfuncs_free();
  }
}

#else
# define pmalloc(s,func,file,line) malloc(s)
# define pcalloc(n,s,func,file,line) calloc(n,s)
# define prealloc(p,s,func,file,line) realloc(p,s)
# define pfree(p,func,file,line) free(p)
# define pstrdup(s,func,file,line) strdup(s)
#endif



typedef struct {
  void * p;
  long max, len;
} Array;

Array * array_new()
{
  Array * a = pmalloc (sizeof(Array),__func__,__FILE__,__LINE__);
  a->p = NULL;
  a->max = a->len = 0;
  return a;
}

void array_free (Array * a)
{
  if (a->max > 0)
    pfree (a->p,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
}

void array_append (Array * a, void * elem, size_t size)
{
  if (a->len + size >= a->max) {
    a->max += max (size, 4096);
    a->p = prealloc (a->p, a->max,__func__,__FILE__,__LINE__);
  }
  memcpy (((char *)a->p) + a->len, elem, size);
  a->len += size;
}



#if TRACE == 1
#include <extrae_user_events.h>

typedef struct {

  Array index, stack;
  extrae_type_t type;
} Trace;

Trace trace_func = {
  {NULL, sizeof(char *), 0, 0}, {NULL, sizeof(int), 0, 0},
  60000010,
};

Trace trace_mpi_func = {
  {NULL, sizeof(char *), 0, 0}, {NULL, sizeof(int), 0, 0},
  60000011,
};

static int lookup_func (Array * a, const char * func)
{
  for (int i = 0; i < a->len; i++) {
    char * s = ((char **)a->p)[i];
    if (!strcmp (func, s))
      return i + 1;
  }
  char * s = pstrdup (func,__func__,__FILE__,__LINE__);
  array_append (a, &s);
  return a->len;
}

static void trace_push (Trace * t, const char * func)
{
  int value = lookup_func (&t->index, func);
  Extrae_eventandcounters (t->type, value);
  array_append (&t->stack, &value);
}

static void trace_pop (Trace * t, const char * func)
{
  assert (t->stack.len > 0);
  t->stack.len--;
  int value = t->stack.len > 0 ? ((int *)t->stack.p)[t->stack.len - 1] : 0;
  Extrae_eventandcounters (t->type, value);
}

static void trace_define (Trace * t, char * description)
{
  if (t->index.len > 0) {
    extrae_value_t values[t->index.len + 1];
    char * names[t->index.len + 1], ** func = (char **) t->index.p;
    names[0] = "OTHER";
    values[0] = 0;
    unsigned len = 1;
    for (int i = 0; i < t->index.len; i++, func++) {
      names[len] = *func;
      values[len++] = i + 1;
    }
    Extrae_define_event_type (&t->type, description, &len, values, names);
  }
}

static void trace_free (Trace * t)
{
  char ** func = (char **) t->index.p;
  for (int i = 0; i < t->index.len; i++, func++)
    pfree (*func,__func__,__FILE__,__LINE__);
  pfree (t->index.p,__func__,__FILE__,__LINE__);
  pfree (t->stack.p,__func__,__FILE__,__LINE__);
}

static void trace_off()
{
  trace_define (&trace_func, "Basilisk functions");
  trace_define (&trace_mpi_func, "Basilisk functions (MPI-related)");
  trace_free (&trace_func);
  trace_free (&trace_mpi_func);
}






# define trace(func, file, line) trace_push (&trace_func, func)
# define end_trace(func, file, line) trace_pop (&trace_func, func)

#elif TRACE

typedef struct {
  char * func, * file;
  int line, calls;
  double total, self;
} TraceIndex;

struct {
  Array stack, index;
  double t0;
} Trace = {
  {NULL, 0, 0}, {NULL, 0, 0},
  -1
};

static void trace_add (const char * func, const char * file, int line,
         double total, double self)
{
  TraceIndex * t = (TraceIndex *) Trace.index.p;
  int i, len = Trace.index.len/sizeof(TraceIndex);
  for (i = 0; i < len; i++, t++)
    if (t->line == line && !strcmp (func, t->func) && !strcmp (file, t->file))
      break;
  if (i == len) {
    TraceIndex t = {pstrdup(func,__func__,__FILE__,__LINE__), pstrdup(file,__func__,__FILE__,__LINE__), line, 1, total, self};
    array_append (&Trace.index, &t, sizeof(TraceIndex));
  }
  else
    t->calls++, t->total += total, t->self += self;
}

static void trace (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  if (Trace.t0 < 0)
    Trace.t0 = tv.tv_sec + tv.tv_usec/1e6;
  double t[2] = {(tv.tv_sec - Trace.t0) + tv.tv_usec/1e6, 0.};
  array_append (&Trace.stack, t, 2*sizeof(double));




}

static void end_trace (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  double te = (tv.tv_sec - Trace.t0) + tv.tv_usec/1e6;
  double * t = (double *) Trace.stack.p;
  assert (Trace.stack.len >= 2*sizeof(double));
  t += Trace.stack.len/sizeof(double) - 2;
  Trace.stack.len -= 2*sizeof(double);
  double dt = te - t[0];




  trace_add (func, file, line, dt, dt - t[1]);
  if (Trace.stack.len >= 2*sizeof(double)) {
    t -= 2;
    t[1] += dt;
  }
}

static int compar_self (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  return t1->self < t2->self;
}

static void trace_off()
{
  int i, len = Trace.index.len/sizeof(TraceIndex);
  double total = 0.;
  TraceIndex * t;
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    total += t->self;
  qsort (Trace.index.p, len, sizeof(TraceIndex), compar_self);
  fprintf (qstderr(), "   calls    total     self   %% total   function\n");
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++) {
    fprintf (qstderr(), "%8d   %6.2f   %6.2f     %4.1f%%   %s():%s:%d\n",
      t->calls, t->total, t->self, t->self*100./total,
      t->func, t->file, t->line);
    pfree (t->func,__func__,__FILE__,__LINE__); pfree (t->file,__func__,__FILE__,__LINE__);
  }

  pfree (Trace.index.p,__func__,__FILE__,__LINE__);
  Trace.index.p = NULL;
  Trace.index.len = Trace.index.max = 0;

  pfree (Trace.stack.p,__func__,__FILE__,__LINE__);
  Trace.stack.p = NULL;
  Trace.stack.len = Trace.stack.max = 0;
}

#else
# define trace(...)
# define end_trace(...)
#endif



#if _OPENMP

#include <omp.h>
#define OMP(x) _Pragma(#x)
#define tid() omp_get_thread_num()
#define pid() 0
#define npe() omp_get_num_threads()
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_double(v,op)

#elif _MPI

#include <mpi.h>
#define OMP(x)

static bool in_prof = false;
static double prof_start, _prof;
#define prof_start(name)\
  assert (!in_prof); in_prof = true;\
  prof_start = MPI_Wtime();\

#line 557

#define prof_stop()\
  assert (in_prof); in_prof = false;\
  _prof = MPI_Wtime();\
  mpi_time += _prof - prof_start;\

#line 562


#if FAKE_MPI
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_double(v,op)
#else

static int mpi_all_reduce0 (void *sendbuf, void *recvbuf, int count,
       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{ trace ("mpi_all_reduce0", "/home/popinet/basilisk/src/common.h", 571);
  { int _ret =  MPI_Allreduce (sendbuf, recvbuf, count, datatype, op, comm); end_trace("mpi_all_reduce0", "/home/popinet/basilisk/src/common.h", 572);  return _ret; }
 end_trace("mpi_all_reduce0", "/home/popinet/basilisk/src/common.h", 573); }
#define mpi_all_reduce(v,type,op) {\
  prof_start ("mpi_all_reduce");\
  union { int a; float b; double c;} global;\
  mpi_all_reduce0 (&(v), &global, 1, type, op, MPI_COMM_WORLD);\
  memcpy (&(v), &global, sizeof (v));\
  prof_stop();\
}\

#line 581

#define mpi_all_reduce_double(v,op) {\
  prof_start ("mpi_all_reduce");\
  double global, tmp = v;\
  mpi_all_reduce0 (&tmp, &global, 1, MPI_DOUBLE, op, MPI_COMM_WORLD);\
  v = global;\
  prof_stop();\
}\

#line 589


#endif

static int mpi_rank, mpi_npe;
#define tid() mpi_rank
#define pid() mpi_rank
#define npe() mpi_npe

#define QFILE FILE

static FILE * qstderr (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "log-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systderr;
  }
  return fp;
}

static FILE * qstdout (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "out-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systdout;
  }
  return fp;
}

static void finalize (void)
{
  MPI_Finalize();
}

void mpi_init()
{
  int initialized;
  MPI_Initialized (&initialized);
  if (!initialized) {
    MPI_Init (NULL, NULL);
    MPI_Comm_set_errhandler (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
    atexit (finalize);
    MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &mpi_npe);
    srand (mpi_rank + 1);
    if (mpi_rank > 0) {
      ferr = fopen ("/dev/null", "w");
      fout = fopen ("/dev/null", "w");
    }
    else {
      ferr = qstderr();
      fout = qstdout();
    }
    char * etrace = getenv ("MALLOC_TRACE"), name[80];
    if (etrace && mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      setenv ("MALLOC_TRACE", name, 1);
    }
#if MTRACE == 1
    etrace = getenv ("MTRACE");
    if (!etrace)
      etrace = "mtrace";
    if (mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      pmtrace.fp = fopen (name, "w");
      pmtrace.fname = systrdup(name);
    }
    else {
      pmtrace.fp = fopen (etrace, "w");
      pmtrace.fname = systrdup(etrace);
    }
#endif
  }
}

#else

#define OMP(x)
#define tid() 0
#define pid() 0
#define npe() 1
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_double(v,op)

#endif

void init_solver()
{
#if _MPI
  mpi_init();
#elif MTRACE == 1
  char * etrace = getenv ("MTRACE");
  pmtrace.fp = fopen (etrace ? etrace : "mtrace", "w");
  pmtrace.fname = systrdup (etrace ? etrace : "mtrace");
#endif
}



#define OMP_PARALLEL() OMP(omp parallel) { _OMPSTART
#define OMP_END_PARALLEL() _OMPEND }
#define _OMPSTART
#define _OMPEND

#define NOT_UNUSED(x) (void)(x)

#define VARIABLES ;
#define val(a,k,l,m) data(k,l,m)[a.i]

double _val_higher_dimension = 0.;
#define _val_higher_dimension(x,a,b,c) _val_higher_dimension
#line 720 "/home/popinet/basilisk/src/common.h"
#if (_GNU_SOURCE || __APPLE__)
double undefined;
# if __APPLE__
# include <stdint.h>
# include "fp_osx.h"
# endif
# define enable_fpe(flags) feenableexcept (flags)
# define disable_fpe(flags) fedisableexcept (flags)
static void set_fpe (void) {
  int64_t lnan = 0x7ff0000000000001;
  assert (sizeof (int64_t) == sizeof (double));
  memcpy (&undefined, &lnan, sizeof (double));
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
#else
# define undefined DBL_MAX
# define enable_fpe(flags)
# define disable_fpe(flags)
static void set_fpe (void) {}
#endif


void * grid = NULL;

double X0 = 0., Y0 = 0., Z0 = 0.;

double L0 = 1.;




int N = 16;


typedef struct { int i; } scalar;

typedef struct {
  scalar x;

  scalar y;


  scalar z;

} vector;

typedef struct {
  vector x;

  vector y;


  vector z;

} tensor;

typedef struct {
  double x, y, z;
} coord;
#line 791 "/home/popinet/basilisk/src/common.h"
struct _origin { double x, y, z; };

void origin (struct _origin p) {
  X0 = p.x; Y0 = p.y; Z0 = p.z;
}

void size (double L) {
  L0 = L;
}

double zero (double s0, double s1, double s2) { return 0.; }
#line 810 "/home/popinet/basilisk/src/common.h"
  enum { right, left, top, bottom, front, back };

int nboundary = 2*3;



#define dirichlet(x) (2.*(x) - val(_s,0,0,0))
#define dirichlet_homogeneous() (- val(_s,0,0,0))
#define neumann(x) (Delta*(x) + val(_s,0,0,0))
#define neumann_homogeneous() (val(_s,0,0,0))

double * _constant = NULL;
extern size_t datasize;
typedef struct _Point Point;

#line 1 "/home/popinet/basilisk/src/grid/boundaries.h"


typedef struct _Boundary Boundary;

struct _Boundary {
  void (* destroy) (Boundary * b);
  void (* level) (const Boundary * b, scalar * list, int l);

  void (* restriction) (const Boundary * b, scalar * list, int l);

  void (* halo_restriction) (const Boundary * b, scalar * list, int l);
  void (* halo_prolongation) (const Boundary * b, scalar * list,
      int l, int depth);
};

void no_halo_restriction (const Boundary * b, scalar * list, int l) {}

static Boundary ** boundaries = NULL;

void add_boundary (Boundary * b) {
  int len = 0;
  if (boundaries) {
    Boundary ** i = boundaries;
    while (*i++) len++;
  }
  boundaries = prealloc (boundaries, sizeof(Boundary *)*(len + 2),__func__,__FILE__,__LINE__);
  boundaries[len] = b;
  boundaries[len+1] = NULL;
}

void free_boundaries() {
  if (!boundaries)
    return;
  Boundary ** i = boundaries, * b;
  while ((b = *i++))
    if (b->destroy)
      b->destroy (b);
    else
      pfree (b,__func__,__FILE__,__LINE__);
  pfree (boundaries,__func__,__FILE__,__LINE__);
  boundaries = NULL;
}
#line 53 "/home/popinet/basilisk/src/grid/boundaries.h"
typedef struct {
  Boundary parent;
  int d;
} BoxBoundary;
#line 826 "/home/popinet/basilisk/src/common.h"



typedef struct {

#line 831 "/home/popinet/basilisk/src/common.h"

  double (** boundary) (Point, Point, scalar);
  double (** boundary_homogeneous) (Point, Point, scalar);
  double (* gradient) (double, double, double);
  void (* delete) (scalar);
  char * name;
  struct {
    int x;

    int y;


    int z;

  } d;
  vector v;
  bool face;

#line 17 "/home/popinet/basilisk-octree/src/grid/multigrid-common.h"

  void (* prolongation) (Point, scalar);
  void (* coarsen) (Point, scalar);

#line 7 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"

  void (* refine) (Point, scalar);

#line 76 "/home/popinet/basilisk-octree/src/fractions.h"

  vector n;

#line 14 "/home/popinet/basilisk-octree/src/tension.h"

  double sigma;
  scalar kappa;

} _Attributes;
_Attributes * _attribute;
#line 829 "/home/popinet/basilisk/src/common.h"























int list_len (scalar * list)
{
  if (!list) return 0;
  int ns = 0;
  if (list) for (scalar s = *list, *_i0 = list; ((scalar *)&s)->i >= 0; s = *++_i0) ns++;
  return ns;
}

scalar * list_append (scalar * list, scalar s)
{
  int len = list_len (list);
  list = prealloc (list, sizeof (scalar)*(len + 2),__func__,__FILE__,__LINE__);
  list[len] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_add (scalar * list, scalar s)
{
  if (list) for (scalar t = *list, *_i1 = list; ((scalar *)&t)->i >= 0; t = *++_i1)
    if (t.i == s.i)
      return list;
  return list_append (list, s);
}

int list_lookup (scalar * l, scalar s)
{
  if (l != NULL)
    if (l) for (scalar s1 = *l, *_i2 = l; ((scalar *)&s1)->i >= 0; s1 = *++_i2)
      if (s1.i == s.i)
 return true;
  return false;
}

scalar * list_copy (scalar * l)
{
  scalar * list = NULL;
  if (l != NULL)
    if (l) for (scalar s = *l, *_i3 = l; ((scalar *)&s)->i >= 0; s = *++_i3)
      list = list_append (list, s);
  return list;
}

scalar * list_concat (scalar * l1, scalar * l2)
{
  scalar * l3 = list_copy (l1);
  if (l2) for (scalar s = *l2, *_i4 = l2; ((scalar *)&s)->i >= 0; s = *++_i4)
    l3 = list_append (l3, s);
  return l3;
}

void list_print (scalar * l, FILE * fp)
{
  int i = 0;
  if (l) for (scalar s = *l, *_i5 = l; ((scalar *)&s)->i >= 0; s = *++_i5)
    fprintf (fp, "%s%s", i++ == 0 ? "{" : ",", _attribute[s.i].name);
  fputs (i > 0 ? "}\n" : "{}\n", fp);
}

int vectors_len (vector * list)
{
  if (!list) return 0;
  int nv = 0;
  if (list) for (vector v = *list, *_i6 = list; ((scalar *)&v)->i >= 0; v = *++_i6) nv++;
  return nv;
}

vector * vectors_append (vector * list, vector v)
{
  int len = vectors_len (list);
  list = prealloc (list, sizeof (vector)*(len + 2),__func__,__FILE__,__LINE__);
  list[len] = v;
  list[len + 1] = (vector){{-1}};
  return list;
}

vector * vectors_add (vector * list, vector v)
{
  if (list) for (vector w = *list, *_i7 = list; ((scalar *)&w)->i >= 0; w = *++_i7) {
    bool id = true;
    {
#line 932

      if (w.x.i != v.x.i)
 id = false;
#line 932

      if (w.y.i != v.y.i)
 id = false;
#line 932

      if (w.z.i != v.z.i)
 id = false;}
    if (id)
      return list;
  }
  return vectors_append (list, v);
}

vector * vectors_copy (vector * l)
{
  vector * list = NULL;
  if (l != NULL)
    if (l) for (vector v = *l, *_i8 = l; ((scalar *)&v)->i >= 0; v = *++_i8)
      list = vectors_append (list, v);
  return list;
}

vector * vectors_from_scalars (scalar * s)
{
  vector * list = NULL;
  while (s->i >= 0) {
    vector v;
    {
#line 955
 {
      assert (s->i >= 0);
      v.x = *s++;
    }
#line 955
 {
      assert (s->i >= 0);
      v.y = *s++;
    }
#line 955
 {
      assert (s->i >= 0);
      v.z = *s++;
    }}
    list = vectors_append (list, v);
  }
  return list;
}

int tensors_len (tensor * list)
{
  if (!list) return 0;
  int nt = 0;
  if (list) for (tensor t = *list, *_i9 = list; ((scalar *)&t)->i >= 0; t = *++_i9) nt++;
  return nt;
}

tensor * tensors_append (tensor * list, tensor t)
{
  int len = tensors_len (list);
  list = prealloc (list, sizeof (tensor)*(len + 2),__func__,__FILE__,__LINE__);
  list[len] = t;
  list[len + 1] = (tensor){{{-1}}};
  return list;
}

tensor * tensors_from_vectors (vector * v)
{
  tensor * list = NULL;
  while (v->x.i >= 0) {
    tensor t;
    {
#line 986
 {
      assert (v->x.i >= 0);
      t.x = *v++;
    }
#line 986
 {
      assert (v->y.i >= 0);
      t.y = *v++;
    }
#line 986
 {
      assert (v->z.i >= 0);
      t.z = *v++;
    }}
    list = tensors_append (list, t);
  }
  return list;
}

scalar * all = NULL;



scalar (* init_scalar) (scalar, const char *);
vector (* init_vector) (vector, const char *);
tensor (* init_tensor) (tensor, const char *);
vector (* init_face_vector) (vector, const char *);





typedef struct _Event Event;
typedef int (* Expr) (int *, double *, Event *);

struct _Event {
  int last, nexpr;
  int (* action) (const int, const double, Event *);
  Expr expr[3];
  int * arrayi;
  double * arrayt;
  char * file;
  int line;
  char * name;
  double t;
  int i, a;
  void * data;
};

static Event * Events = NULL;

double tnext = 0;
void init_events (void);
void event_register (Event event);
void _init_solver (void);



#if _MPI
static double mpi_time = 0.;
#endif

typedef struct {
  clock_t c;
  struct timeval tv;
  double tm;
} timer;

timer timer_start (void)
{
  timer t;
  t.c = clock();
  gettimeofday (&t.tv, NULL);
#if _MPI
  t.tm = mpi_time;
#endif
  return t;
}

double timer_elapsed (timer t)
{
  struct timeval tvend;
  gettimeofday (&tvend, NULL);
  return ((tvend.tv_sec - t.tv.tv_sec) +
   (tvend.tv_usec - t.tv.tv_usec)/1e6);
}



vector zerof= {{_NVARMAX + 0},{_NVARMAX + 1},{_NVARMAX + 2}};
vector unityf= {{_NVARMAX + 3},{_NVARMAX + 4},{_NVARMAX + 5}};
scalar unity= {_NVARMAX + 6};
scalar zeroc= {_NVARMAX + 7};



 vector fm = {{_NVARMAX + 3},{_NVARMAX + 4},{_NVARMAX + 5}};
 scalar cm = {(_NVARMAX + 6)};



static FILE ** qpopen_pipes = NULL;

FILE * qpopen (const char * command, const char * type)
{
  FILE * fp = popen (command, type);
  if (fp) {
    FILE ** i = qpopen_pipes;
    int n = 0;
    while (i && *i) { n++; i++; }
    qpopen_pipes = prealloc (qpopen_pipes, sizeof(FILE *)*(n + 2),__func__,__FILE__,__LINE__);
    qpopen_pipes[n] = fp;
    qpopen_pipes[n+1] = NULL;
  }
  return fp;
}

int qpclose (FILE * fp)
{
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i == fp)
      *i = (FILE *) 1;
    i++;
  }
  return pclose (fp);
}

static void qpclose_all()
{
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i != (FILE *) 1)
      pclose (*i);
    i++;
  }
  pfree (qpopen_pipes,__func__,__FILE__,__LINE__);
  qpopen_pipes = NULL;
}






FILE * lfopen (const char * name, const char * mode)
{
  char fname[80];
  sprintf (fname, "%s-%d", name, pid());
  return fopen (fname, mode);
}



void * matrix_new (int n, int p, size_t size)
{
  void ** m = pmalloc (n*sizeof (void *),__func__,__FILE__,__LINE__);
  char * a = pmalloc (n*p*size,__func__,__FILE__,__LINE__);
  for (int i = 0; i < n; i++)
    m[i] = a + i*p*size;
  return m;
}

double matrix_inverse (double ** m, int n, double pivmin)
{
  int indxc[n], indxr[n], ipiv[n];
  int i, icol = 0, irow = 0, j, k, l, ll;
  double big, dum, pivinv, minpiv = HUGE;

  for (j = 0; j < n; j++)
    ipiv[j] = -1;

  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if (ipiv[j] != 0)
 for (k = 0; k < n; k++) {
   if (ipiv[k] == -1) {
     if (fabs (m[j][k]) >= big) {
       big = fabs (m[j][k]);
       irow = j;
       icol = k;
     }
   }
 }
    ipiv[icol]++;
    if (irow != icol)
      for (l = 0; l < n; l++)
 swap (double, m[irow][l], m[icol][l]);
    indxr[i] = irow;
    indxc[i] = icol;
    if (fabs (m[icol][icol]) <= pivmin)
      return 0.;
    if (fabs (m[icol][icol]) < minpiv)
      minpiv = fabs (m[icol][icol]);
    pivinv = 1.0/m[icol][icol];
    m[icol][icol] = 1.0;
    for (l = 0; l < n; l++) m[icol][l] *= pivinv;
    for (ll = 0; ll < n; ll++)
      if (ll != icol) {
 dum = m[ll][icol];
 m[ll][icol] = 0.0;
 for (l = 0; l < n; l++)
   m[ll][l] -= m[icol][l]*dum;
      }
  }
  for (l = n - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
 swap (double, m[k][indxr[l]], m[k][indxc[l]]);
  }
  return minpiv;
}

void matrix_free (void * m)
{
  pfree (((void **) m)[0],__func__,__FILE__,__LINE__);
  pfree (m,__func__,__FILE__,__LINE__);
}
#line 14 "atomisation-cpp.c"
#line 1 "grid/octree.h"
#line 1 "/home/popinet/basilisk-octree/src/grid/octree.h"


#line 1 "grid/tree.h"
#line 1 "/home/popinet/basilisk-octree/src/grid/tree.h"
#line 1 "grid/mempool.h"
#line 1 "/home/popinet/basilisk-octree/src/grid/mempool.h"





typedef struct _Pool Pool;

struct _Pool {
  Pool * next;
};

typedef struct {
  char * first, * lastb;
  size_t size;
  size_t poolsize;
  Pool * pool, * last;
} Mempool;

typedef struct {
  char * next;
} FreeBlock;

Mempool * mempool_new (size_t poolsize, size_t size)
{

  assert (poolsize % 8 == 0);
  assert (size >= sizeof(FreeBlock));


  poolsize = min(1 << 20, poolsize + sizeof(Pool));
  Mempool * m = pcalloc (1, sizeof(Mempool),__func__,__FILE__,__LINE__);
  m->poolsize = poolsize;
  m->size = size;
  return m;
}

void mempool_destroy (Mempool * m)
{
  Pool * p = m->pool;
  while (p) {
    Pool * next = p->next;
    pfree (p,__func__,__FILE__,__LINE__);
    p = next;
  }
  pfree (m,__func__,__FILE__,__LINE__);
}

void * mempool_alloc (Mempool * m)
{
  if (!m->first) {

    Pool * p = pmalloc (m->poolsize,__func__,__FILE__,__LINE__);
    p->next = NULL;
    if (m->last)
      m->last->next = p;
    else
      m->pool = p;
    m->last = p;
    m->first = m->lastb = ((char *)m->last) + sizeof(Pool);
    FreeBlock * b = (FreeBlock *) m->first;
    b->next = NULL;
  }
  void * ret = m->first;
  FreeBlock * b = ret;
  char * next = b->next;
  if (!next) {
    m->lastb += m->size;
    next = m->lastb;
    if (next + m->size > ((char *) m->last) + m->poolsize)
      next = NULL;
    else {
      FreeBlock * b = (FreeBlock *) next;
      b->next = NULL;
    }
  }
  m->first = next;
#if TRASH
  double * v = ret;
  for (int i = 0; i < m->size/sizeof(double); i++)
    v[i] = undefined;
#endif
  return ret;
}

void * mempool_alloc0 (Mempool * m)
{
  void * ret = mempool_alloc (m);
  memset (ret, 0, m->size);
  return ret;
}

void mempool_free (Mempool * m, void * p)
{
#if TRASH
  double * v = p;
  for (int i = 0; i < m->size/sizeof(double); i++)
    v[i] = undefined;
#endif
  FreeBlock * b = p;
  b->next = m->first;
  m->first = p;
}
#line 2 "/home/popinet/basilisk-octree/src/grid/tree.h"
#line 22 "/home/popinet/basilisk-octree/src/grid/tree.h"
typedef struct {
  unsigned short flags;
  unsigned short neighbors;
  int pid;
} Cell;

enum {
  active = 1 << 0,
  leaf = 1 << 1,
  halo = 1 << 2,
  border = 1 << 3,
  vertex = 1 << 4,
  user = 5,

  face_x = 1 << 0

  , face_y = 1 << 1


  , face_z = 1 << 2

};

#define is_active(cell) ((cell).flags & active)
#define is_leaf(cell) ((cell).flags & leaf)
#define is_coarse() ((cell).neighbors > 0)
#define is_border(cell) ((cell).flags & border)
#define is_local(cell) ((cell).pid == pid())



typedef struct {
  int i;

  int j;


  int k;

} IndexLevel;

typedef struct {
  IndexLevel * p;
  int n, nm;
} CacheLevel;

typedef struct {
  int i;

  int j;


  int k;

  int level, flags;
} Index;

typedef struct {
  Index * p;
  int n, nm;
} Cache;




static void * new_refarray (size_t len, size_t size) {
  return pcalloc (len + 1, size,__func__,__FILE__,__LINE__);
}

static void refarray (void * p, size_t len, size_t size) {
  int * refcount = (int *)(((char *)p) + len*size);
  (*refcount)++;
}

static bool unrefarray (void * p, size_t len, size_t size) {
  int * refcount = (int *)(((char *)p) + len*size);
  (*refcount)--;
  if (*refcount == 0) {
    pfree (p,__func__,__FILE__,__LINE__);
    return true;
  }
  return false;
}




typedef struct {





  char **** m;

  Mempool * pool;
  int nc;
  int len;
} Layer;

static size_t _size (size_t depth)
{
  return (1 << depth) + 2*2;
}

static size_t poolsize (size_t depth, size_t size)
{






  return cube(_size(depth))*size;

}

static Layer * new_layer (int depth)
{
  Layer * l = pmalloc (sizeof (Layer),__func__,__FILE__,__LINE__);
  l->len = _size (depth);
  if (depth == 0)
    l->pool = NULL;
  else {
    size_t size = sizeof(Cell) + datasize;


    l->pool = mempool_new (poolsize (depth, size), (1 << 3)*size);
  }
  l->m = pcalloc (l->len, sizeof (char *),__func__,__FILE__,__LINE__);
  l->nc = 0;
  return l;
}

static void destroy_layer (Layer * l)
{
  if (l->pool)
    mempool_destroy (l->pool);
  pfree (l->m,__func__,__FILE__,__LINE__);
  pfree (l,__func__,__FILE__,__LINE__);
}
#line 180 "/home/popinet/basilisk-octree/src/grid/tree.h"
static void layer_add_row (Layer * l, int i, int j)
{
  if (!l->m[i]) {
    l->m[i] = new_refarray (l->len, sizeof (char *));
    l->nc++;
  }
  refarray (l->m[i], l->len, sizeof(char *));

  if (!l->m[i][j])
    l->m[i][j] = new_refarray (l->len, sizeof (char *));
  refarray (l->m[i][j], l->len, sizeof(char *));

}

static bool layer_remove_row (Layer * l, int i, int j)
{

  if (unrefarray (l->m[i][j], l->len, sizeof (char *)))
    l->m[i][j] = NULL;

  if (unrefarray (l->m[i], l->len, sizeof (char *))) {
    l->m[i] = NULL;
    if (--l->nc == 0) {
      destroy_layer (l);
      return true;
    }
    assert (l->nc >= 0);
  }
  return false;
}




typedef struct {
  Layer ** L;
  int depth;

  Cache leaves;
  Cache faces;
  Cache vertices;
  Cache refined;
  CacheLevel * active;
  CacheLevel * prolongation;
  CacheLevel * restriction;
  CacheLevel * boundary;

  bool dirty;
} Quadtree;



struct _Point {

  int i;

  int j;


  int k;

  int level;
};
static Point last_point;



static void cache_level_append (CacheLevel * c, Point p)
{
  if (c->n >= c->nm) {
    c->nm += 128;
    c->p = prealloc (c->p, sizeof (IndexLevel)*c->nm,__func__,__FILE__,__LINE__);
  }
  c->p[c->n].i = p.i;

  c->p[c->n].j = p.j;


  c->p[c->n].k = p.k;

  c->n++;
}

static void cache_level_shrink (CacheLevel * c)
{
  if (c->nm > (c->n/128 + 1)*128) {
    c->nm = (c->n/128 + 1)*128;
    assert (c->nm > c->n);
    c->p = prealloc (c->p, sizeof (Index)*c->nm,__func__,__FILE__,__LINE__);
  }
}

static void cache_append (Cache * c, Point p, unsigned short flags)
{
  if (c->n >= c->nm) {
    c->nm += 128;
    c->p = prealloc (c->p, sizeof (Index)*c->nm,__func__,__FILE__,__LINE__);
  }
  c->p[c->n].i = p.i;

  c->p[c->n].j = p.j;


  c->p[c->n].k = p.k;

  c->p[c->n].level = p.level;
  c->p[c->n].flags = flags;
  c->n++;
}

void cache_shrink (Cache * c)
{
  cache_level_shrink ((CacheLevel *)c);
}
#line 349 "/home/popinet/basilisk-octree/src/grid/tree.h"
#define allocated(a,l,n) (point.i+a >= 0 &&\
         point.i+a < (1 << point.level) + 2*2 &&\
         ((Quadtree *)grid)->L[point.level]->m[point.i+a] &&\
         point.j+l >= 0 &&\
         point.j+l < (1 << point.level) + 2*2 &&\
         ((Quadtree *)grid)->L[point.level]->m[point.i+a][point.j+l] &&\
         point.k+n >= 0 &&\
         point.k+n < (1 << point.level) + 2*2 &&\
         ((Quadtree *)grid)->L[point.level]->m[point.i+a][point.j+l]\
         [point.k+n])\

#line 359


#define NEIGHBOR(a,l,n) (((Quadtree *)grid)->L[point.level]->m[point.i+a][point.j+l]\
                                          [point.k+n])\

#line 363

#define PARENT(a,l,n) (((Quadtree *)grid)->L[point.level-1]->m[(point.i+2)/2+a]\
    [(point.j+2)/2+l][(point.k+2)/2+n])\

#line 366

#define allocated_child(a,l,n) (level < depth() &&\
         point.i > 0 && point.i <= (1 << level) + 2 &&\
         point.j > 0 && point.j <= (1 << level) + 2 &&\
         point.k > 0 && point.k <= (1 << level) + 2 &&\
         ((Quadtree *)grid)->L[point.level+1]->m[2*point.i-2 +a]\
   && ((Quadtree *)grid)->L[point.level+1]->m[2*point.i-2 +a][2*point.j-2 +l]\
   && ((Quadtree *)grid)->L[point.level+1]->m[2*point.i-2 +a][2*point.j-2 +l]\
         [2*point.k-2 +n])\

#line 375

#define CHILD(a,l,n) (((Quadtree *)grid)->L[point.level+1]->m[2*point.i-2 +a]\
      [2*point.j-2 +l][2*point.k-2 +n])\

#line 378


#define CELL(m) (*((Cell *)(m)))


#define depth() (((Quadtree *)grid)->depth)
#define aparent(k,l,n) CELL(PARENT(k,l,n))
#define child(k,l,n) CELL(CHILD(k,l,n))


#define cell CELL(NEIGHBOR(0,0,0))
#define neighbor(k,l,n) CELL(NEIGHBOR(k,l,n))
#define neighborp(l,m,n) (Point) {\
    point.i + l,\
\
    point.j + m,\
\
\
    point.k + n,\
\
    point.level }\

#line 399



#define data(k,l,n) ((double *) (NEIGHBOR(k,l,n) + sizeof(Cell)))
#define fine(a,k,l,n) ((double *) (CHILD(k,l,n) + sizeof(Cell)))[a.i]
#define coarse(a,k,l,n) ((double *) (PARENT(k,l,n) + sizeof(Cell)))[a.i]

#define POINT_VARIABLES\
  VARIABLES\
  int level = point.level; NOT_UNUSED(level);\
\
\
\
\
\
\
\
  struct { int x, y, z; } child = {\
    2*((point.i+2)%2)-1, 2*((point.j+2)%2)-1, 2*((point.k+2)%2)-1\
  };\
\
  NOT_UNUSED(child);\
  Point parent = point; NOT_UNUSED(parent);\
  parent.level--;\
  parent.i = (point.i + 2)/2;\
\
  parent.j = (point.j + 2)/2;\
\
\
  parent.k = (point.k + 2)/2;\
\

#line 430


#line 1 "grid/foreach_cell.h"
#line 1 "/home/popinet/basilisk-octree/src/grid/foreach_cell.h"
#line 66 "/home/popinet/basilisk-octree/src/grid/foreach_cell.h"
#define foreach_cell_root(root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point;\
\
\
\
\
\
    int kg = 0; NOT_UNUSED(kg);\
    struct { int l, i, j, k, stage; } stack[20];\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].k = root.k; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; point.k = stack[_s].k; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
 POINT_VARIABLES;\
\

#line 89

#define end_foreach_cell_root()\
        if (point.level < ((Quadtree *)grid)->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 1; };\
          { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
        }\
        break;\
      }\
\
      case 1: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 2; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
      case 2: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 3; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;\
      case 3: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 4; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
      case 4: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 5; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;\
      case 5: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 6; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
      case 6: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 7; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;\
      case 7: { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
\
      }\
    }\
  }\

#line 123


#define foreach_cell() {\
\
\
\
\
\
  Point root = {2,2,2,0};\
\
  foreach_cell_root (root)\

#line 134

#define end_foreach_cell() end_foreach_cell_root() }

#define foreach_cell_all() {\
  Point root = { .level = 0 };\
  for (root.i = 0; root.i <= 2*2; root.i++)\
\
    for (root.j = 0; root.j <= 2*2; root.j++)\
\
\
      for (root.k = 0; root.k <= 2*2; root.k++)\
\
 foreach_cell_root (root)\

#line 147

#define end_foreach_cell_all() end_foreach_cell_root() }

#define foreach_cell_post_root(condition, root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point;\
\
\
\
\
\
    int kg = 0; NOT_UNUSED(kg);\
    struct { int l, i, j, k, stage; } stack[20];\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].k = root.k; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; point.k = stack[_s].k; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
        POINT_VARIABLES;\
 if (point.level == ((Quadtree *)grid)->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 8; };\
 }\
 else {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 1; };\
   if (condition)\
     { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 }\
 break;\
      }\
\
      case 1:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 2; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 2:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 3; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 break;\
      case 3:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 4; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 4:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 5; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 break;\
      case 5:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 6; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 6:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 7; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 break;\
      case 7:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 8; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
 break;\
\
      default: {\
        POINT_VARIABLES;\
\

#line 244

#define end_foreach_cell_post_root()\
      }\
      }\
    }\
  }\

#line 250


#define foreach_cell_post(condition)\
  {\
\
\
\
\
\
    Point root = {2,2,2,0};\
\
    foreach_cell_post_root(condition, root)\

#line 262

#define end_foreach_cell_post() end_foreach_cell_post_root() }

#define foreach_cell_post_all(condition) {\
  Point root = { .level = 0 };\
  for (root.i = 0; root.i <= 2*2; root.i++)\
\
    for (root.j = 0; root.j <= 2*2; root.j++)\
\
\
      for (root.k = 0; root.k <= 2*2; root.k++)\
\
 foreach_cell_post_root (condition, root)\

#line 275

#define end_foreach_cell_post_all() end_foreach_cell_post_root() }
#line 433 "/home/popinet/basilisk-octree/src/grid/tree.h"
#line 468 "/home/popinet/basilisk-octree/src/grid/tree.h"
#define foreach_child() {\
  int _i = 2*point.i - 2, _j = 2*point.j - 2, _k = 2*point.k - 2;\
  point.level++;\
  for (int _l = 0; _l < 2; _l++) {\
    point.i = _i + _l;\
    for (int _m = 0; _m < 2; _m++) {\
      point.j = _j + _m;\
      for (int _n = 0; _n < 2; _n++) {\
 point.k = _k + _n;\
 POINT_VARIABLES;\

#line 478

#define end_foreach_child()\
      }\
    }\
  }\
  point.i = (_i + 2)/2;point.j = (_j + 2)/2;point.k = (_k + 2)/2;\
  point.level--;\
}\

#line 486

#define foreach_child_break() _l = _m = _n = 2
#line 496 "/home/popinet/basilisk-octree/src/grid/tree.h"
#define foreach_cache(_cache,clause) {\
  int ig = 0; NOT_UNUSED(ig);\
\
  int jg = 0; NOT_UNUSED(jg);\
\
\
  int kg = 0; NOT_UNUSED(kg);\
\
  OMP_PARALLEL()\
\
\
\
\
\
  Point point = {2,2,2,0};\
\
  int _k; unsigned short _flags; NOT_UNUSED(_flags);\
  OMP(omp for schedule(static) clause)\
  for (_k = 0; _k < _cache.n; _k++) {\
    point.i = _cache.p[_k].i;\
\
    point.j = _cache.p[_k].j;\
\
\
    point.k = _cache.p[_k].k;\
\
    point.level = _cache.p[_k].level;\
    _flags = _cache.p[_k].flags;\
    POINT_VARIABLES;\

#line 525

#define end_foreach_cache() } OMP_END_PARALLEL() }

#define foreach_cache_level(_cache,_l,clause) {\
  int ig = 0; NOT_UNUSED(ig);\
\
  int jg = 0; NOT_UNUSED(jg);\
\
\
  int kg = 0; NOT_UNUSED(kg);\
\
  OMP_PARALLEL()\
\
\
\
\
\
  Point point = {2,2,2,0};\
\
  point.level = _l;\
  int _k;\
  OMP(omp for schedule(static) clause)\
  for (_k = 0; _k < _cache.n; _k++) {\
    point.i = _cache.p[_k].i;\
\
    point.j = _cache.p[_k].j;\
\
\
    point.k = _cache.p[_k].k;\
\
    POINT_VARIABLES;\

#line 556

#define end_foreach_cache_level() } OMP_END_PARALLEL() }

#define foreach_boundary(_l) {\
  if (_l <= depth()) {\
    { if (((Quadtree *)grid)->dirty) update_cache_f(); };\
    CacheLevel _boundary = ((Quadtree *)grid)->boundary[_l];\
    foreach_cache_level (_boundary,_l,)\

#line 564

#define end_foreach_boundary() end_foreach_cache_level(); }}

#define foreach_halo(_name,_l) {\
  if (_l <= depth()) {\
    { if (((Quadtree *)grid)->dirty) update_cache_f(); };\
    CacheLevel _cache = ((Quadtree *)grid)->_name[_l];\
    foreach_cache_level (_cache, _l,)\

#line 572

#define end_foreach_halo() end_foreach_cache_level(); }}

#line 1 "grid/neighbors.h"
#line 1 "/home/popinet/basilisk-octree/src/grid/neighbors.h"
#line 35 "/home/popinet/basilisk-octree/src/grid/neighbors.h"
#define foreach_neighbor(_s) {\
  int _nn = _s + 0 ? _s + 0 : 2;\
  int _i = point.i, _j = point.j, _k = point.k;\
  for (int _l = - _nn; _l <= _nn; _l++) {\
    point.i = _i + _l;\
    for (int _m = - _nn; _m <= _nn; _m++) {\
      point.j = _j + _m;\
      for (int _n = - _nn; _n <= _nn; _n++) {\
 point.k = _k + _n;\
 POINT_VARIABLES;\

#line 45

#define end_foreach_neighbor()\
      }\
    }\
  }\
  point.i = _i; point.j = _j; point.k = _k;\
}\

#line 52

#define foreach_neighbor_break() _l = _m = _n = _nn + 1
#line 576 "/home/popinet/basilisk-octree/src/grid/tree.h"

static inline bool has_local_children (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 578 "/home/popinet/basilisk-octree/src/grid/tree.h"

   { foreach_child()
    if (is_local(cell))
      return true; end_foreach_child(); }
  return false;
}

static inline void cache_append_face (Point point, unsigned short flags)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 586 "/home/popinet/basilisk-octree/src/grid/tree.h"

  Quadtree * q = grid;
  cache_append (&q->faces, point, flags);
#line 600 "/home/popinet/basilisk-octree/src/grid/tree.h"
  {
#line 600

    if (flags & face_x)
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)
   if (!(neighbor(0,i,j).flags & vertex)) {
     cache_append (&q->vertices, neighborp(0,i,j), 0);
     neighbor(0,i,j).flags |= vertex;
   }
#line 600

    if (flags & face_y)
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)
   if (!(neighbor(j,0,i).flags & vertex)) {
     cache_append (&q->vertices, neighborp(j,0,i), 0);
     neighbor(j,0,i).flags |= vertex;
   }
#line 600

    if (flags & face_z)
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)
   if (!(neighbor(i,j,0).flags & vertex)) {
     cache_append (&q->vertices, neighborp(i,j,0), 0);
     neighbor(i,j,0).flags |= vertex;
   }}

}

static void update_cache_f (void)
{
  Quadtree * q = grid;

   { foreach_cache(q->vertices,){

#line 615 "/home/popinet/basilisk-octree/src/grid/tree.h"

    if (level <= depth() && allocated(0,0,0))
      cell.flags &= ~vertex; } end_foreach_cache(); }


  q->leaves.n = q->faces.n = q->vertices.n = 0;
  for (int l = 0; l <= depth(); l++)
    q->active[l].n = q->prolongation[l].n =
      q->restriction[l].n = q->boundary[l].n = 0;

  const unsigned short fboundary = 1 << user;
   { foreach_cell(){

#line 626 "/home/popinet/basilisk-octree/src/grid/tree.h"
 {
    if (is_local(cell)) {

      assert (is_active(cell));
      cache_level_append (&q->active[level], point);
    }

    if (!(cell.pid < 0))

       { foreach_neighbor (2)
 if (allocated(0,0,0) && (cell.pid < 0) && !(cell.flags & fboundary)) {
   cache_level_append (&q->boundary[level], point);
   cell.flags |= fboundary;
 } end_foreach_neighbor(); }
    if (is_leaf (cell)) {
      if (is_local(cell)) {
 cache_append (&q->leaves, point, 0);

 unsigned short flags = 0;
 {
#line 645

   if ((neighbor(-1,0,0).pid < 0) || (!is_leaf(neighbor(-1,0,0)) && !neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) ||
       is_leaf(neighbor(-1,0,0)))
     flags |= face_x;
#line 645

   if ((neighbor(0,-1,0).pid < 0) || (!is_leaf(neighbor(0,-1,0)) && !neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) ||
       is_leaf(neighbor(0,-1,0)))
     flags |= face_y;
#line 645

   if ((neighbor(0,0,-1).pid < 0) || (!is_leaf(neighbor(0,0,-1)) && !neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0) ||
       is_leaf(neighbor(0,0,-1)))
     flags |= face_z;}
 if (flags)
   cache_append (&q->faces, point, flags);
 {
#line 651

   if ((neighbor(1,0,0).pid < 0) || (!is_leaf(neighbor(1,0,0)) && !neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0) ||
       (!is_local(neighbor(1,0,0)) && is_leaf(neighbor(1,0,0))))
     cache_append (&q->faces, neighborp(1,0,0), face_x);
#line 651

   if ((neighbor(0,1,0).pid < 0) || (!is_leaf(neighbor(0,1,0)) && !neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0) ||
       (!is_local(neighbor(0,1,0)) && is_leaf(neighbor(0,1,0))))
     cache_append (&q->faces, neighborp(0,1,0), face_y);
#line 651

   if ((neighbor(0,0,1).pid < 0) || (!is_leaf(neighbor(0,0,1)) && !neighbor(0,0,1).neighbors && neighbor(0,0,1).pid >= 0) ||
       (!is_local(neighbor(0,0,1)) && is_leaf(neighbor(0,0,1))))
     cache_append (&q->faces, neighborp(0,0,1), face_z);}

 for (int i = 0; i <= 1; i++)

   for (int j = 0; j <= 1; j++)


     for (int k = 0; k <= 1; k++)

       if (!(neighbor(i,j,k).flags & vertex)) {
  cache_append (&q->vertices, neighborp(i,j,k), 0);
  neighbor(i,j,k).flags |= vertex;
       }

        if (cell.neighbors > 0)
   cache_level_append (&q->prolongation[level], point);
 cell.flags &= ~halo;
      }
      else if (!(cell.pid < 0) || is_local(aparent(0,0,0))) {

 unsigned short flags = 0;
 {
#line 675

   if (allocated(-1,0,0) &&
       is_local(neighbor(-1,0,0)) && (!is_leaf(neighbor(-1,0,0)) && !neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0))
     flags |= face_x;
#line 675

   if (allocated(0,-1,0) &&
       is_local(neighbor(0,-1,0)) && (!is_leaf(neighbor(0,-1,0)) && !neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0))
     flags |= face_y;
#line 675

   if (allocated(0,0,-1) &&
       is_local(neighbor(0,0,-1)) && (!is_leaf(neighbor(0,0,-1)) && !neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0))
     flags |= face_z;}
 if (flags)
   cache_append_face (point, flags);
 {
#line 681

   if (allocated(1,0,0) && is_local(neighbor(1,0,0)) &&
       (!is_leaf(neighbor(1,0,0)) && !neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0))
     cache_append_face (neighborp(1,0,0), face_x);
#line 681

   if (allocated(0,1,0) && is_local(neighbor(0,1,0)) &&
       (!is_leaf(neighbor(0,1,0)) && !neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0))
     cache_append_face (neighborp(0,1,0), face_y);
#line 681

   if (allocated(0,0,1) && is_local(neighbor(0,0,1)) &&
       (!is_leaf(neighbor(0,0,1)) && !neighbor(0,0,1).neighbors && neighbor(0,0,1).pid >= 0))
     cache_append_face (neighborp(0,0,1), face_z);}
      }
      continue;
    }
    else {
#if _MPI

      if (is_local(cell) && (cell.flags & halo))
 cache_level_append (&q->restriction[level], point);
#else
      bool restriction = level > 0 && (aparent(0,0,0).flags & halo);
      if (!restriction) {

 if (is_active(cell)) {
    { foreach_neighbor()
     if (allocated(0,0,0) && is_leaf(cell) && !(cell.pid < 0))
       restriction = true, foreach_neighbor_break(); end_foreach_neighbor(); }
 }
 else
    { foreach_neighbor()
     if (allocated(0,0,0) && is_leaf(cell) && is_local(cell))
       restriction = true, foreach_neighbor_break(); end_foreach_neighbor(); }
      }
      if (restriction) {

 cell.flags |= halo;
 if (is_local(cell))
   cache_level_append (&q->restriction[level], point);
      }
      else
 cell.flags &= ~halo;
#endif
    }
  } } end_foreach_cell(); }


  cache_shrink (&q->leaves);
  cache_shrink (&q->faces);
  cache_shrink (&q->vertices);
  for (int l = 0; l <= depth(); l++) {
    cache_level_shrink (&q->active[l]);
    cache_level_shrink (&q->prolongation[l]);
    cache_level_shrink (&q->restriction[l]);
    cache_level_shrink (&q->boundary[l]);
  }

  q->dirty = false;

  for (int l = depth(); l >= 0; l--) {
     { foreach_boundary (l){

#line 733 "/home/popinet/basilisk-octree/src/grid/tree.h"

      cell.flags &= ~(fboundary|halo); } end_foreach_boundary(); }

     { foreach_halo (restriction, l){

#line 736 "/home/popinet/basilisk-octree/src/grid/tree.h"

       { foreach_child()
        if ((cell.pid < 0))
   cell.flags |= halo; end_foreach_child(); }; } end_foreach_halo(); }
  }
}

#define foreach(clause) { if (((Quadtree *)grid)->dirty) update_cache_f(); }; foreach_cache(((Quadtree *)grid)->leaves, clause)
#define end_foreach() end_foreach_cache()

#define foreach_face_generic(clause)\
  { if (((Quadtree *)grid)->dirty) update_cache_f(); };\
  foreach_cache(((Quadtree *)grid)->faces, clause) 
#line 747

#define end_foreach_face_generic() end_foreach_cache()

#define is_face_x() (_flags & face_x)

#define is_face_y() (_flags & face_y)


#define is_face_z() (_flags & face_z)


#define foreach_vertex(clause)\
  { if (((Quadtree *)grid)->dirty) update_cache_f(); };\
  foreach_cache(((Quadtree *)grid)->vertices, clause) {\
    x -= Delta/2.;\
\
    y -= Delta/2.;\
\
\
    z -= Delta/2.;\
\

#line 769

#define end_foreach_vertex() } end_foreach_cache()
#line 781 "/home/popinet/basilisk-octree/src/grid/tree.h"
#define foreach_fine_to_coarse(clause) {\
  { if (((Quadtree *)grid)->dirty) update_cache_f(); };\
  for (int _l = depth() - 1; _l >= 0; _l--) {\
    CacheLevel _active = ((Quadtree *)grid)->active[_l];\
    foreach_cache_level (_active,_l,clause)\
      if (!is_leaf (cell)) {\

#line 787

#define end_foreach_fine_to_coarse() } end_foreach_cache_level(); } }

#define foreach_level(l) {\
  if (l <= depth()) {\
    { if (((Quadtree *)grid)->dirty) update_cache_f(); };\
    CacheLevel _active = ((Quadtree *)grid)->active[l];\
    foreach_cache_level (_active,l,)\

#line 795

#define end_foreach_level() end_foreach_cache_level(); }}

#define foreach_coarse_level(l) foreach_level(l) if (!is_leaf(cell)) {
#define end_foreach_coarse_level() } end_foreach_level()

#define foreach_level_or_leaf(l) {\
  for (int _l1 = l; _l1 >= 0; _l1--)\
    foreach_level(_l1)\
      if (_l1 == l || is_leaf (cell)) {\

#line 805

#define end_foreach_level_or_leaf() } end_foreach_level(); }

#define foreach_leaf() foreach_cell()\
  if (is_leaf (cell)) {\
    if (is_active(cell) && is_local(cell)) {\

#line 811

#define end_foreach_leaf() } continue; } end_foreach_cell()

#if TRASH
# undef trash
# define trash quadtree_trash
#endif

void quadtree_trash (void * alist)
{
  scalar * list = alist;
  Quadtree * q = grid;

  for (int l = 0; l <= q->depth; l++) {
    Layer * L = q->L[l];
    for (int i = 0; i < L->len; i++)
      if (L->m[i])





 for (int j = 0; j < L->len; j++)
   if (L->m[i][j])





          for (int k = 0; k < L->len; k++)
     if (L->m[i][j][k])
       if (list) for (scalar s = *list, *_i10 = list; ((scalar *)&s)->i >= 0; s = *++_i10)
           if (!is_constant(s))
    ((double *)(L->m[i][j][k] + sizeof(Cell)))[s.i] = undefined;


  }
}

#define cache_level_resize(name, a)\
{\
  for (int i = 0; i <= q->depth - a; i++)\
    pfree (q->name[i].p,__func__,__FILE__,__LINE__);\
  pfree (q->name,__func__,__FILE__,__LINE__);\
  q->name = pcalloc (q->depth + 1, sizeof (CacheLevel),__func__,__FILE__,__LINE__);\
}\

#line 857


static void update_depth (int inc)
{
  Quadtree * q = ((Quadtree *)grid);
  q->depth += inc;
  q->L = &(q->L[-1]);
  q->L = prealloc(q->L, sizeof (Layer *)*(q->depth + 2),__func__,__FILE__,__LINE__);
  q->L = &(q->L[1]);
  if (inc > 0)
    q->L[q->depth] = new_layer (q->depth);
  cache_level_resize (active, inc);
  cache_level_resize (prolongation, inc);
  cache_level_resize (restriction, inc);
  cache_level_resize (boundary, inc);
}

static void alloc_children (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 875 "/home/popinet/basilisk-octree/src/grid/tree.h"

  if (point.level == ((Quadtree *)grid)->depth)
    update_depth (+1);
  else if (allocated_child(0,0,0))
    return;
#line 907 "/home/popinet/basilisk-octree/src/grid/tree.h"
  Layer * L = ((Quadtree *)grid)->L[point.level + 1];
  size_t len = sizeof(Cell) + datasize;
  char * b = mempool_alloc0 (L->pool);
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++) {
      layer_add_row (L, 2*point.i - 2 + k, 2*point.j - 2 + l);
      for (int n = 0; n < 2; n++) {
 assert (!CHILD(k,l,n));
 CHILD(k,l,n) = b;
 b += len;
      }
    }


  int pid = cell.pid;
   { foreach_child() {
    cell.pid = pid;
#if TRASH
    if (all) for (scalar s = *all, *_i11 = all; ((scalar *)&s)->i >= 0; s = *++_i11)
      val(s,0,0,0) = undefined;
#endif
  } end_foreach_child(); }
}
#line 963 "/home/popinet/basilisk-octree/src/grid/tree.h"
static void free_children (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 964 "/home/popinet/basilisk-octree/src/grid/tree.h"


  Layer * L = ((Quadtree *)grid)->L[point.level + 1];
  assert (CHILD(0,0,0));
  mempool_free (L->pool, CHILD(0,0,0));
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++) {
      for (int n = 0; n < 2; n++)
 CHILD(k,l,n) = NULL;
      if (layer_remove_row (L, 2*point.i - 2 + k, 2*point.j - 2 + l)) {
 assert (point.level + 1 == ((Quadtree *)grid)->depth);
 update_depth (-1);
      }
    }
}


void increment_neighbors (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 982 "/home/popinet/basilisk-octree/src/grid/tree.h"

  ((Quadtree *)grid)->dirty = true;
   { foreach_neighbor (2/2)
    if (cell.neighbors++ == 0)
      alloc_children (point); end_foreach_neighbor(); }
}

void decrement_neighbors (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 990 "/home/popinet/basilisk-octree/src/grid/tree.h"

  ((Quadtree *)grid)->dirty = true;
   { foreach_neighbor (2/2)
    if (allocated(0,0,0)) {
      cell.neighbors--;
      if (cell.neighbors == 0)
 free_children (point);
    } end_foreach_neighbor(); }
  if (cell.neighbors) {
    int pid = cell.pid;
     { foreach_child() {
      cell.flags = 0;
      cell.pid = pid;
    } end_foreach_child(); }
  }
}

void realloc_scalar (void)
{

  Quadtree * q = grid;
  size_t newlen = sizeof(Cell) + datasize;
  size_t oldlen = newlen - sizeof(double);

  size_t len = _size(0);
  for (int i = 0; i < len; i++)



    for (int j = 0; j < len; j++)



      for (int k = 0; k < len; k++)
 q->L[0]->m[i][j][k] = prealloc (q->L[0]->m[i][j][k], newlen,__func__,__FILE__,__LINE__);



  for (int l = 1; l <= q->depth; l++) {
    Layer * L = q->L[l];
    size_t len = L->len;
    Mempool * oldpool = L->pool;
    L->pool = mempool_new (poolsize (l, newlen), (1 << 3)*newlen);
    for (int i = 0; i < len; i += 2)
      if (L->m[i]) {
#line 1043 "/home/popinet/basilisk-octree/src/grid/tree.h"
 for (int j = 0; j < len; j += 2)
   if (L->m[i][j]) {
#line 1054 "/home/popinet/basilisk-octree/src/grid/tree.h"
     for (int k = 0; k < len; k += 2)
       if (L->m[i][j][k]) {
  char * new = mempool_alloc (L->pool);
  for (int l = 0; l < 2; l++)
    for (int m = 0; m < 2; m++)
      for (int n = 0; n < 2; n++) {
        memcpy (new, L->m[i+l][j+m][k+n], oldlen);
        L->m[i+l][j+m][k+n] = new;
        new += newlen;
      }
       }

   }

      }
    mempool_destroy (oldpool);
  }
}





#define VN v.x
#define VT v.y





#if _MPI
# define disable_fpe_for_mpi() disable_fpe (FE_DIVBYZERO|FE_INVALID)
# define enable_fpe_for_mpi() enable_fpe (FE_DIVBYZERO|FE_INVALID)
#else
# define disable_fpe_for_mpi()
# define enable_fpe_for_mpi()
#endif

static inline void no_coarsen (Point point, scalar s);

static bool normal_neighbor (Point point, bool (cond)(Point),
        scalar * scalars, vector * vectors)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1096 "/home/popinet/basilisk-octree/src/grid/tree.h"

  for (int k = 1; k <= 2; k++)
    {
#line 1098

      for (int i = -k; i <= k; i += 2*k)
 if ((allocated(i,0,0) && !(neighbor(i,0,0).pid < 0) && cond(neighborp(i,0,0)))) {
   Point neighbor = neighborp(i,0,0);
   int id = (- cell.pid - 1);
   if (scalars) for (scalar s = *scalars, *_i12 = scalars; ((scalar *)&s)->i >= 0; s = *++_i12)
     val(s,0,0,0) = _attribute[s.i].boundary[id](neighbor, point, s);
   if (vectors) for (vector v = *vectors, *_i13 = vectors; ((scalar *)&v)->i >= 0; v = *++_i13) {
     scalar vn = VN;
     val(v.x,0,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.x);

     scalar vt = VT;
     val(v.y,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.y);


     val(v.z,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.z);

   }
   return true;
 }
#line 1098

      for (int i = -k; i <= k; i += 2*k)
 if ((allocated(0,i,0) && !(neighbor(0,i,0).pid < 0) && cond(neighborp(0,i,0)))) {
   Point neighbor = neighborp(0,i,0);
   int id = (- cell.pid - 1);
   if (scalars) for (scalar s = *scalars, *_i12 = scalars; ((scalar *)&s)->i >= 0; s = *++_i12)
     val(s,0,0,0) = _attribute[s.i].boundary[id](neighbor, point, s);
   if (vectors) for (vector v = *vectors, *_i13 = vectors; ((scalar *)&v)->i >= 0; v = *++_i13) {
     scalar vn = VN;
     val(v.y,0,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.y);

     scalar vt = VT;
     val(v.z,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.z);


     val(v.x,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.x);

   }
   return true;
 }
#line 1098

      for (int i = -k; i <= k; i += 2*k)
 if ((allocated(0,0,i) && !(neighbor(0,0,i).pid < 0) && cond(neighborp(0,0,i)))) {
   Point neighbor = neighborp(0,0,i);
   int id = (- cell.pid - 1);
   if (scalars) for (scalar s = *scalars, *_i12 = scalars; ((scalar *)&s)->i >= 0; s = *++_i12)
     val(s,0,0,0) = _attribute[s.i].boundary[id](neighbor, point, s);
   if (vectors) for (vector v = *vectors, *_i13 = vectors; ((scalar *)&v)->i >= 0; v = *++_i13) {
     scalar vn = VN;
     val(v.z,0,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.z);

     scalar vt = VT;
     val(v.x,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.x);


     val(v.y,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.y);

   }
   return true;
 }}
  return false;
}

static bool diagonal_neighbor_2D (Point point, bool (cond)(Point),
      scalar * scalars, vector * vectors)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1123 "/home/popinet/basilisk-octree/src/grid/tree.h"


  for (int k = 1; k <= 2; k++)

    {
#line 1127


      for (int i = -k; i <= k; i += 2*k)
 for (int j = -k; j <= k; j += 2*k)
   if (allocated(i,j,0) && (allocated(i,j,0) && !(neighbor(i,j,0).pid < 0) && cond(neighborp(i,j,0))) &&
       allocated(i,0,0) && (neighbor(i,0,0).pid < 0) &&
       allocated(0,j,0) && (neighbor(0,j,0).pid < 0)) {
     Point n = neighborp(i,j,0),
       n1 = neighborp(i,0,0), n2 = neighborp(0,j,0);
     int id1 = (- neighbor(i,0,0).pid - 1), id2 = (- neighbor(0,j,0).pid - 1);
     if (scalars) for (scalar s = *scalars, *_i14 = scalars; ((scalar *)&s)->i >= 0; s = *++_i14)
       val(s,0,0,0) = (_attribute[s.i].boundary[id1](n,n1,s) + _attribute[s.i].boundary[id2](n,n2,s) -
       val(s,i,j,0));
     if (vectors) for (vector v = *vectors, *_i15 = vectors; ((scalar *)&v)->i >= 0; v = *++_i15) {
       scalar vt = VT, vn = VN;
       val(v.x,0,0,0) = (_attribute[vt.i].boundary[id1](n,n1,v.x) +
         _attribute[vn.i].boundary[id2](n,n2,v.x) -
         val(v.x,i,j,0));
       val(v.y,0,0,0) = (_attribute[vn.i].boundary[id1](n,n1,v.y) +
         _attribute[vt.i].boundary[id2](n,n2,v.y) -
         val(v.y,i,j,0));

       val(v.z,0,0,0) = (_attribute[vt.i].boundary[id1](n,n1,v.z) +
         _attribute[vt.i].boundary[id2](n,n2,v.z) -
         val(v.z,i,j,0));

     }
     return true;
   }
#line 1127


      for (int i = -k; i <= k; i += 2*k)
 for (int j = -k; j <= k; j += 2*k)
   if (allocated(0,i,j) && (allocated(0,i,j) && !(neighbor(0,i,j).pid < 0) && cond(neighborp(0,i,j))) &&
       allocated(0,i,0) && (neighbor(0,i,0).pid < 0) &&
       allocated(0,0,j) && (neighbor(0,0,j).pid < 0)) {
     Point n = neighborp(0,i,j),
       n1 = neighborp(0,i,0), n2 = neighborp(0,0,j);
     int id1 = (- neighbor(0,i,0).pid - 1), id2 = (- neighbor(0,0,j).pid - 1);
     if (scalars) for (scalar s = *scalars, *_i14 = scalars; ((scalar *)&s)->i >= 0; s = *++_i14)
       val(s,0,0,0) = (_attribute[s.i].boundary[id1](n,n1,s) + _attribute[s.i].boundary[id2](n,n2,s) -
       val(s,0,i,j));
     if (vectors) for (vector v = *vectors, *_i15 = vectors; ((scalar *)&v)->i >= 0; v = *++_i15) {
       scalar vt = VT, vn = VN;
       val(v.y,0,0,0) = (_attribute[vt.i].boundary[id1](n,n1,v.y) +
         _attribute[vn.i].boundary[id2](n,n2,v.y) -
         val(v.y,0,i,j));
       val(v.z,0,0,0) = (_attribute[vn.i].boundary[id1](n,n1,v.z) +
         _attribute[vt.i].boundary[id2](n,n2,v.z) -
         val(v.z,0,i,j));

       val(v.x,0,0,0) = (_attribute[vt.i].boundary[id1](n,n1,v.x) +
         _attribute[vt.i].boundary[id2](n,n2,v.x) -
         val(v.x,0,i,j));

     }
     return true;
   }
#line 1127


      for (int i = -k; i <= k; i += 2*k)
 for (int j = -k; j <= k; j += 2*k)
   if (allocated(j,0,i) && (allocated(j,0,i) && !(neighbor(j,0,i).pid < 0) && cond(neighborp(j,0,i))) &&
       allocated(0,0,i) && (neighbor(0,0,i).pid < 0) &&
       allocated(j,0,0) && (neighbor(j,0,0).pid < 0)) {
     Point n = neighborp(j,0,i),
       n1 = neighborp(0,0,i), n2 = neighborp(j,0,0);
     int id1 = (- neighbor(0,0,i).pid - 1), id2 = (- neighbor(j,0,0).pid - 1);
     if (scalars) for (scalar s = *scalars, *_i14 = scalars; ((scalar *)&s)->i >= 0; s = *++_i14)
       val(s,0,0,0) = (_attribute[s.i].boundary[id1](n,n1,s) + _attribute[s.i].boundary[id2](n,n2,s) -
       val(s,j,0,i));
     if (vectors) for (vector v = *vectors, *_i15 = vectors; ((scalar *)&v)->i >= 0; v = *++_i15) {
       scalar vt = VT, vn = VN;
       val(v.z,0,0,0) = (_attribute[vt.i].boundary[id1](n,n1,v.z) +
         _attribute[vn.i].boundary[id2](n,n2,v.z) -
         val(v.z,j,0,i));
       val(v.x,0,0,0) = (_attribute[vn.i].boundary[id1](n,n1,v.x) +
         _attribute[vt.i].boundary[id2](n,n2,v.x) -
         val(v.x,j,0,i));

       val(v.y,0,0,0) = (_attribute[vt.i].boundary[id1](n,n1,v.y) +
         _attribute[vt.i].boundary[id2](n,n2,v.y) -
         val(v.y,j,0,i));

     }
     return true;
   }}

  return false;
}

static bool diagonal_neighbor_3D (Point point, bool (cond)(Point),
      scalar * scalars, vector * vectors)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1162 "/home/popinet/basilisk-octree/src/grid/tree.h"


  for (int n = 1; n <= 2; n++)
    for (int i = -n; i <= n; i += 2*n)
      for (int j = -n; j <= n; j += 2*n)
 for (int k = -n; k <= n; k += 2*n)
   if ((allocated(i,j,k) && !(neighbor(i,j,k).pid < 0) && cond(neighborp(i,j,k))) &&
       (neighbor(i,j,0).pid < 0) &&
       (neighbor(i,0,k).pid < 0) &&
       (neighbor(0,j,k).pid < 0)) {
     Point
       n0 = neighborp(i,j,k),
       n1 = neighborp(i,j,0),
       n2 = neighborp(i,0,k),
       n3 = neighborp(0,j,k);
     int
       id1 = (- neighbor(i,j,0).pid - 1),
       id2 = (- neighbor(i,0,k).pid - 1),
       id3 = (- neighbor(0,j,k).pid - 1);
     if (scalars) for (scalar s = *scalars, *_i16 = scalars; ((scalar *)&s)->i >= 0; s = *++_i16)
       val(s,0,0,0) = (_attribute[s.i].boundary[id1](n0,n1,s) +
       _attribute[s.i].boundary[id2](n0,n2,s) +
       _attribute[s.i].boundary[id3](n0,n3,s) -
       2.*val(s,i,j,k));
     if (vectors) for (vector v = *vectors, *_i17 = vectors; ((scalar *)&v)->i >= 0; v = *++_i17) {
       scalar vt = VT, vn = VN;
       val(v.x,0,0,0) = (_attribute[vt.i].boundary[id1](n0,n1,v.x) +
         _attribute[vt.i].boundary[id2](n0,n2,v.x) +
         _attribute[vn.i].boundary[id3](n0,n3,v.x) -
         2.*val(v.x,i,j,k));
       val(v.y,0,0,0) = (_attribute[vt.i].boundary[id1](n0,n1,v.y) +
         _attribute[vn.i].boundary[id2](n0,n2,v.y) +
         _attribute[vt.i].boundary[id3](n0,n3,v.y) -
         2.*val(v.y,i,j,k));
       val(v.z,0,0,0) = (_attribute[vn.i].boundary[id1](n0,n1,v.z) +
         _attribute[vt.i].boundary[id2](n0,n2,v.z) +
         _attribute[vt.i].boundary[id3](n0,n3,v.z) -
         2.*val(v.z,i,j,k));
     }
     return true;
   }

  return false;
}



#line 1208

static Point tangential_neighbor_x (Point point, bool (cond)(Point))
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1210 "/home/popinet/basilisk-octree/src/grid/tree.h"

  for (int k = 1; k <= 2; k++)
    for (int j = -k; j <= k; j += 2*k) {
      if ((allocated(0,j,0) && !(neighbor(0,j,0).pid < 0) && cond(neighborp(0,j,0))) || (allocated(-1,j,0) && !(neighbor(-1,j,0).pid < 0) && cond(neighborp(-1,j,0))))
 return neighborp(0,j,0);


      if ((allocated(0,0,j) && !(neighbor(0,0,j).pid < 0) && cond(neighborp(0,0,j))) || (allocated(-1,0,j) && !(neighbor(-1,0,j).pid < 0) && cond(neighborp(-1,0,j))))
 return neighborp(0,0,j);

    }
  return (Point){.level = -1};
}
#line 1208

static Point tangential_neighbor_y (Point point, bool (cond)(Point))
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1210 "/home/popinet/basilisk-octree/src/grid/tree.h"

  for (int k = 1; k <= 2; k++)
    for (int j = -k; j <= k; j += 2*k) {
      if ((allocated(0,0,j) && !(neighbor(0,0,j).pid < 0) && cond(neighborp(0,0,j))) || (allocated(0,-1,j) && !(neighbor(0,-1,j).pid < 0) && cond(neighborp(0,-1,j))))
 return neighborp(0,0,j);


      if ((allocated(j,0,0) && !(neighbor(j,0,0).pid < 0) && cond(neighborp(j,0,0))) || (allocated(j,-1,0) && !(neighbor(j,-1,0).pid < 0) && cond(neighborp(j,-1,0))))
 return neighborp(j,0,0);

    }
  return (Point){.level = -1};
}
#line 1208

static Point tangential_neighbor_z (Point point, bool (cond)(Point))
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1210 "/home/popinet/basilisk-octree/src/grid/tree.h"

  for (int k = 1; k <= 2; k++)
    for (int j = -k; j <= k; j += 2*k) {
      if ((allocated(j,0,0) && !(neighbor(j,0,0).pid < 0) && cond(neighborp(j,0,0))) || (allocated(j,0,-1) && !(neighbor(j,0,-1).pid < 0) && cond(neighborp(j,0,-1))))
 return neighborp(j,0,0);


      if ((allocated(0,j,0) && !(neighbor(0,j,0).pid < 0) && cond(neighborp(0,j,0))) || (allocated(0,j,-1) && !(neighbor(0,j,-1).pid < 0) && cond(neighborp(0,j,-1))))
 return neighborp(0,j,0);

    }
  return (Point){.level = -1};
}


void box_boundaries (int l,
       bool (cond1)(Point), bool (cond)(Point),
       scalar * list)
{
  scalar * scalars = NULL;
  vector * vectors = NULL, * faces = NULL;
  if (list) for (scalar s = *list, *_i18 = list; ((scalar *)&s)->i >= 0; s = *++_i18)
    if (!is_constant(s) && _attribute[s.i].refine != no_coarsen) {
      if (_attribute[s.i].v.x.i == s.i) {
 if (_attribute[s.i].face)
   faces = vectors_add (faces, _attribute[s.i].v);
 else
   vectors = vectors_add (vectors, _attribute[s.i].v);
      }
      else if (_attribute[s.i].v.x.i < 0)
 scalars = list_add (scalars, s);
    }

   { foreach_boundary (l){

#line 1243 "/home/popinet/basilisk-octree/src/grid/tree.h"

    if (cond1 (point)) {
      if (!normal_neighbor (point, cond, scalars, vectors) &&
   !diagonal_neighbor_2D (point, cond, scalars, vectors) &&
   !diagonal_neighbor_3D (point, cond, scalars, vectors)) {

 if (scalars) for (scalar s = *scalars, *_i19 = scalars; ((scalar *)&s)->i >= 0; s = *++_i19)
   val(s,0,0,0) = undefined;
 if (vectors) for (vector v = *vectors, *_i20 = vectors; ((scalar *)&v)->i >= 0; v = *++_i20)
   {
#line 1252

     val(v.x,0,0,0) = undefined;
#line 1252

     val(v.y,0,0,0) = undefined;
#line 1252

     val(v.z,0,0,0) = undefined;}
      }
      if (faces) {
 int id = (- cell.pid - 1);
 {
#line 1257

   for (int i = -1; i <= 1; i += 2) {

     if ((allocated(i,0,0) && !(neighbor(i,0,0).pid < 0) && cond(neighborp(i,0,0)))) {
       Point neighbor = neighborp(i,0,0);
       if (faces) for (vector v = *faces, *_i21 = faces; ((scalar *)&v)->i >= 0; v = *++_i21) {
  scalar vn = VN;
  if (_attribute[vn.i].boundary[id])
    val(v.x,(i + 1)/2,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.x);
       }
     }

     else if (i == -1) {

       Point neighbor = tangential_neighbor_x (point, cond);
       if (neighbor.level >= 0) {
  if (faces) for (vector v = *faces, *_i22 = faces; ((scalar *)&v)->i >= 0; v = *++_i22) {
    scalar vt = VT;
    val(v.x,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.x);
  }
       }
       else

  if (faces) for (vector v = *faces, *_i23 = faces; ((scalar *)&v)->i >= 0; v = *++_i23)
    val(v.x,0,0,0) = 0.;
     }

   }
#line 1257

   for (int i = -1; i <= 1; i += 2) {

     if ((allocated(0,i,0) && !(neighbor(0,i,0).pid < 0) && cond(neighborp(0,i,0)))) {
       Point neighbor = neighborp(0,i,0);
       if (faces) for (vector v = *faces, *_i21 = faces; ((scalar *)&v)->i >= 0; v = *++_i21) {
  scalar vn = VN;
  if (_attribute[vn.i].boundary[id])
    val(v.y,0,(i + 1)/2,0) = _attribute[vn.i].boundary[id](neighbor, point, v.y);
       }
     }

     else if (i == -1) {

       Point neighbor = tangential_neighbor_y (point, cond);
       if (neighbor.level >= 0) {
  if (faces) for (vector v = *faces, *_i22 = faces; ((scalar *)&v)->i >= 0; v = *++_i22) {
    scalar vt = VT;
    val(v.y,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.y);
  }
       }
       else

  if (faces) for (vector v = *faces, *_i23 = faces; ((scalar *)&v)->i >= 0; v = *++_i23)
    val(v.y,0,0,0) = 0.;
     }

   }
#line 1257

   for (int i = -1; i <= 1; i += 2) {

     if ((allocated(0,0,i) && !(neighbor(0,0,i).pid < 0) && cond(neighborp(0,0,i)))) {
       Point neighbor = neighborp(0,0,i);
       if (faces) for (vector v = *faces, *_i21 = faces; ((scalar *)&v)->i >= 0; v = *++_i21) {
  scalar vn = VN;
  if (_attribute[vn.i].boundary[id])
    val(v.z,0,0,(i + 1)/2) = _attribute[vn.i].boundary[id](neighbor, point, v.z);
       }
     }

     else if (i == -1) {

       Point neighbor = tangential_neighbor_z (point, cond);
       if (neighbor.level >= 0) {
  if (faces) for (vector v = *faces, *_i22 = faces; ((scalar *)&v)->i >= 0; v = *++_i22) {
    scalar vt = VT;
    val(v.z,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.z);
  }
       }
       else

  if (faces) for (vector v = *faces, *_i23 = faces; ((scalar *)&v)->i >= 0; v = *++_i23)
    val(v.z,0,0,0) = 0.;
     }

   }}
      }
    } } end_foreach_boundary(); }

  pfree (scalars,__func__,__FILE__,__LINE__);
  pfree (vectors,__func__,__FILE__,__LINE__);
  pfree (faces,__func__,__FILE__,__LINE__);
}



#undef VN
#undef VT
#define VN _attribute[s.i].v.x
#define VT _attribute[s.i].v.y

static bool retrue (Point point) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1299 "/home/popinet/basilisk-octree/src/grid/tree.h"
 return true; }

static bool retleaf (Point point) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1301 "/home/popinet/basilisk-octree/src/grid/tree.h"
 return is_leaf(cell); }

static bool retleafhalo (Point point) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1304 "/home/popinet/basilisk-octree/src/grid/tree.h"


  return is_leaf(cell) || !cell.neighbors || (cell.flags & halo);
}

static bool retleafhalo1 (Point point) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1309 "/home/popinet/basilisk-octree/src/grid/tree.h"


  return is_leaf(cell) || (cell.flags & halo);
}

static bool rethalo (Point point) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1314 "/home/popinet/basilisk-octree/src/grid/tree.h"


  return cell.flags & halo;
}

static void box_boundary_level (const Boundary * b, scalar * list, int l)
{
  disable_fpe_for_mpi();
  if (l < 0)
    for (l = 0; l <= depth(); l++)
      box_boundaries (l, retrue, retleaf, list);
  else
    box_boundaries (l, retrue, retrue, list);
  enable_fpe_for_mpi();
}

static void box_boundary_halo_restriction (const Boundary * b,
        scalar * list, int l)
{
  disable_fpe_for_mpi();
  if (l == 0)
    return;
  box_boundaries (l, rethalo, retleafhalo1, list);
  enable_fpe_for_mpi();
}

static void box_boundary_halo_prolongation (const Boundary * b,
         scalar * list,
         int l, int depth)
{
  disable_fpe_for_mpi();
#if _MPI
    box_boundaries (l, retrue, retrue, list);
#else
  if (l == depth)
    box_boundaries (l, retrue, retrue, list);
  else
    box_boundaries (l, retrue, retleafhalo, list);
#endif
  enable_fpe_for_mpi();
}
#line 1381 "/home/popinet/basilisk-octree/src/grid/tree.h"
static double periodic_bc (Point, Point, scalar);


#line 1383

static void periodic_boundary_level_x (const Boundary * b, scalar * list, int l)
{
  scalar * list1 = NULL;
  if (list) for (scalar s = *list, *_i24 = list; ((scalar *)&s)->i >= 0; s = *++_i24)
    if (!is_constant(s)) {
      if (_attribute[s.i].face) {

 scalar vt = VT;
 if (_attribute[vt.i].boundary[right] == periodic_bc)
   list1 = list_add (list1, s);

      }
      else if (_attribute[s.i].boundary[right] == periodic_bc) {
 if (_attribute[s.i].v.x.i >= 0)
   {
#line 1398

     list1 = list_add (list1, _attribute[s.i].v.x);
#line 1398

     list1 = list_add (list1, _attribute[s.i].v.y);
#line 1398

     list1 = list_add (list1, _attribute[s.i].v.z);}
 else
   list1 = list_add (list1, s);
      }
    }
  if (!list1)
    return;

  OMP_PARALLEL();



  Point point = {0,0,0};  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1411 "/home/popinet/basilisk-octree/src/grid/tree.h"


  point.level = l < 0 ? depth() : l;
  int n = 1 << point.level;

  int j;
  OMP(omp for schedule(static))
  for (j = 0; j < n + 2*2; j++)


    for (int k = 0; k < n + 2*2; k++)

    {
      for (int i = 0; i < 2; i++)
 if (allocated(i + n,j,k))
   if (list1) for (scalar s = *list1, *_i25 = list1; ((scalar *)&s)->i >= 0; s = *++_i25)
     val(s,i,j,k) = val(s,i + n,j,k);
      for (int i = n + 2; i < n + 2*2; i++)
 if (allocated(i - n,j,k))
   if (list1) for (scalar s = *list1, *_i26 = list1; ((scalar *)&s)->i >= 0; s = *++_i26)
     val(s,i,j,k) = val(s,i - n,j,k);
    }
  OMP_END_PARALLEL();

  pfree (list1,__func__,__FILE__,__LINE__);
}
#line 1383

static void periodic_boundary_level_y (const Boundary * b, scalar * list, int l)
{
  scalar * list1 = NULL;
  if (list) for (scalar s = *list, *_i24 = list; ((scalar *)&s)->i >= 0; s = *++_i24)
    if (!is_constant(s)) {
      if (_attribute[s.i].face) {

 scalar vt = VT;
 if (_attribute[vt.i].boundary[top] == periodic_bc)
   list1 = list_add (list1, s);

      }
      else if (_attribute[s.i].boundary[top] == periodic_bc) {
 if (_attribute[s.i].v.y.i >= 0)
   {
#line 1398

     list1 = list_add (list1, _attribute[s.i].v.y);
#line 1398

     list1 = list_add (list1, _attribute[s.i].v.z);
#line 1398

     list1 = list_add (list1, _attribute[s.i].v.x);}
 else
   list1 = list_add (list1, s);
      }
    }
  if (!list1)
    return;

  OMP_PARALLEL();



  Point point = {0,0,0};  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1411 "/home/popinet/basilisk-octree/src/grid/tree.h"


  point.level = l < 0 ? depth() : l;
  int n = 1 << point.level;

  int j;
  OMP(omp for schedule(static))
  for (j = 0; j < n + 2*2; j++)


    for (int k = 0; k < n + 2*2; k++)

    {
      for (int i = 0; i < 2; i++)
 if (allocated(k,i + n,j))
   if (list1) for (scalar s = *list1, *_i25 = list1; ((scalar *)&s)->i >= 0; s = *++_i25)
     val(s,k,i,j) = val(s,k,i + n,j);
      for (int i = n + 2; i < n + 2*2; i++)
 if (allocated(k,i - n,j))
   if (list1) for (scalar s = *list1, *_i26 = list1; ((scalar *)&s)->i >= 0; s = *++_i26)
     val(s,k,i,j) = val(s,k,i - n,j);
    }
  OMP_END_PARALLEL();

  pfree (list1,__func__,__FILE__,__LINE__);
}
#line 1383

static void periodic_boundary_level_z (const Boundary * b, scalar * list, int l)
{
  scalar * list1 = NULL;
  if (list) for (scalar s = *list, *_i24 = list; ((scalar *)&s)->i >= 0; s = *++_i24)
    if (!is_constant(s)) {
      if (_attribute[s.i].face) {

 scalar vt = VT;
 if (_attribute[vt.i].boundary[front] == periodic_bc)
   list1 = list_add (list1, s);

      }
      else if (_attribute[s.i].boundary[front] == periodic_bc) {
 if (_attribute[s.i].v.z.i >= 0)
   {
#line 1398

     list1 = list_add (list1, _attribute[s.i].v.z);
#line 1398

     list1 = list_add (list1, _attribute[s.i].v.x);
#line 1398

     list1 = list_add (list1, _attribute[s.i].v.y);}
 else
   list1 = list_add (list1, s);
      }
    }
  if (!list1)
    return;

  OMP_PARALLEL();



  Point point = {0,0,0};  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1411 "/home/popinet/basilisk-octree/src/grid/tree.h"


  point.level = l < 0 ? depth() : l;
  int n = 1 << point.level;

  int j;
  OMP(omp for schedule(static))
  for (j = 0; j < n + 2*2; j++)


    for (int k = 0; k < n + 2*2; k++)

    {
      for (int i = 0; i < 2; i++)
 if (allocated(j,k,i + n))
   if (list1) for (scalar s = *list1, *_i25 = list1; ((scalar *)&s)->i >= 0; s = *++_i25)
     val(s,j,k,i) = val(s,j,k,i + n);
      for (int i = n + 2; i < n + 2*2; i++)
 if (allocated(j,k,i - n))
   if (list1) for (scalar s = *list1, *_i26 = list1; ((scalar *)&s)->i >= 0; s = *++_i26)
     val(s,j,k,i) = val(s,j,k,i - n);
    }
  OMP_END_PARALLEL();

  pfree (list1,__func__,__FILE__,__LINE__);
}

#undef VN
#undef VT


#line 1441

static void periodic_boundary_halo_prolongation_x (const Boundary * b,
         scalar * list,
         int l, int depth)
{
  periodic_boundary_level_x (b, list, l);
}
#line 1441

static void periodic_boundary_halo_prolongation_y (const Boundary * b,
         scalar * list,
         int l, int depth)
{
  periodic_boundary_level_y (b, list, l);
}
#line 1441

static void periodic_boundary_halo_prolongation_z (const Boundary * b,
         scalar * list,
         int l, int depth)
{
  periodic_boundary_level_z (b, list, l);
}

static void free_cache (CacheLevel * c)
{
  Quadtree * q = grid;
  for (int l = 0; l <= q->depth; l++)
    pfree (c[l].p,__func__,__FILE__,__LINE__);
  pfree (c,__func__,__FILE__,__LINE__);
}

void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
  Quadtree * q = grid;
  pfree (q->leaves.p,__func__,__FILE__,__LINE__);
  pfree (q->faces.p,__func__,__FILE__,__LINE__);
  pfree (q->vertices.p,__func__,__FILE__,__LINE__);
  pfree (q->refined.p,__func__,__FILE__,__LINE__);


  Layer * L = q->L[0];
#line 1487 "/home/popinet/basilisk-octree/src/grid/tree.h"
  for (int i = 0; i < L->len; i++) {
    for (int j = 0; j < L->len; j++) {
      for (int k = 0; k < L->len; k++)
 pfree (L->m[i][j][k],__func__,__FILE__,__LINE__);
      pfree (L->m[i][j],__func__,__FILE__,__LINE__);
    }
    pfree (L->m[i],__func__,__FILE__,__LINE__);
  }

  for (int l = 1; l <= q->depth; l++) {
    Layer * L = q->L[l];
    for (int i = 0; i < L->len; i++)
      if (L->m[i]) {
 for (int j = 0; j < L->len; j++)
   if (L->m[i][j])
     pfree (L->m[i][j],__func__,__FILE__,__LINE__);
 pfree (L->m[i],__func__,__FILE__,__LINE__);
      }
  }

  for (int l = 0; l <= q->depth; l++)
    destroy_layer (q->L[l]);
  q->L = &(q->L[-1]);
  pfree (q->L,__func__,__FILE__,__LINE__);
  free_cache (q->active);
  free_cache (q->prolongation);
  free_cache (q->restriction);
  free_cache (q->boundary);
  pfree (q,__func__,__FILE__,__LINE__);
  grid = NULL;
}

static void refine_level (int depth);




void init_grid (int n)
{ trace ("init_grid", "/home/popinet/basilisk-octree/src/grid/tree.h", 1525);

  assert (sizeof(Cell) % 8 == 0);

  free_grid();
  int depth = 0;
  while (n > 1) {
    if (n % 2) {
      fprintf (qstderr(), "quadtree: N must be a power-of-two\n");
      exit (1);
    }
    n /= 2;
    depth++;
  }
  Quadtree * q = pcalloc (1, sizeof (Quadtree),__func__,__FILE__,__LINE__);
  q->depth = 0;


  q->L = pmalloc(sizeof (Layer *)*2,__func__,__FILE__,__LINE__);

  q->L[0] = NULL; q->L = &(q->L[1]);

  Layer * L = new_layer (0);
  q->L[0] = L;
#line 1589 "/home/popinet/basilisk-octree/src/grid/tree.h"
  for (int i = 0; i < L->len; i++)
    for (int j = 0; j < L->len; j++) {
      layer_add_row (L, i, j);
      for (int k = 0; k < L->len; k++)
 L->m[i][j][k] = pcalloc (1, sizeof(Cell) + datasize,__func__,__FILE__,__LINE__);
    }
  CELL(L->m[2][2][2]).flags |= leaf;
  if (pid() == 0)
    CELL(L->m[2][2][2]).flags |= active;
  for (int k = -2; k <= 2; k++)
    for (int l = -2; l <= 2; l++)
      for (int n = -2; n <= 2; n++)
 CELL(L->m[2 +k][2 +l][2 +n]).pid =
   (k > 0 ? -1 - right :
    k < 0 ? -1 - left :
    l > 0 ? -1 - top :
    l < 0 ? -1 - bottom :
    n > 0 ? -1 - front :
    n < 0 ? -1 - back :
    0);

  q->active = pcalloc (1, sizeof (CacheLevel),__func__,__FILE__,__LINE__);
  q->prolongation = pcalloc (1, sizeof (CacheLevel),__func__,__FILE__,__LINE__);
  q->restriction = pcalloc (1, sizeof (CacheLevel),__func__,__FILE__,__LINE__);
  q->boundary = pcalloc (1, sizeof (CacheLevel),__func__,__FILE__,__LINE__);
  q->dirty = true;
  grid = q;
  N = 1 << depth;
#if _MPI
  void mpi_boundary_new();
  mpi_boundary_new();
#endif
  { if (((Quadtree *)grid)->dirty) update_cache_f(); };

  Boundary * b = pcalloc (1, sizeof (Boundary),__func__,__FILE__,__LINE__);
  b->level = b->restriction = box_boundary_level;
  b->halo_restriction = box_boundary_halo_restriction;
  b->halo_prolongation = box_boundary_halo_prolongation;
  add_boundary (b);

  {
#line 1629
 {
    Boundary * b = pcalloc (1, sizeof (Boundary),__func__,__FILE__,__LINE__);
    b->level = b->restriction = periodic_boundary_level_x;
    b->halo_prolongation = periodic_boundary_halo_prolongation_x;
    add_boundary (b);
  }
#line 1629
 {
    Boundary * b = pcalloc (1, sizeof (Boundary),__func__,__FILE__,__LINE__);
    b->level = b->restriction = periodic_boundary_level_y;
    b->halo_prolongation = periodic_boundary_halo_prolongation_y;
    add_boundary (b);
  }
#line 1629
 {
    Boundary * b = pcalloc (1, sizeof (Boundary),__func__,__FILE__,__LINE__);
    b->level = b->restriction = periodic_boundary_level_z;
    b->halo_prolongation = periodic_boundary_halo_prolongation_z;
    add_boundary (b);
  }}
  refine_level (depth);
  trash (all);







 end_trace("init_grid", "/home/popinet/basilisk-octree/src/grid/tree.h", 1644); }
#line 1680 "/home/popinet/basilisk-octree/src/grid/tree.h"
struct _locate { double x, y, z; };

Point locate (struct _locate p)
{
  for (int l = depth(); l >= 0; l--) {
    Point point = { .level = l };
    int n = 1 << point.level;
    point.i = (p.x - X0)/L0*n + 2;

    point.j = (p.y - Y0)/L0*n + 2;


    point.k = (p.z - Z0)/L0*n + 2;

    if (point.i >= 0 && point.i < n + 2*2

 && point.j >= 0 && point.j < n + 2*2


 && point.k >= 0 && point.k < n + 2*2

 ) {
      if (allocated(0,0,0) && is_local(cell) && is_leaf(cell))
 return point;
    }
    else
      break;
  }
  Point point = { .level = -1 };
  return point;
}

#line 1 "grid/quadtree-common.h"
#line 1 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"


#line 1 "grid/multigrid-common.h"
#line 1 "/home/popinet/basilisk-octree/src/grid/multigrid-common.h"


#line 1 "grid/cartesian-common.h"
#line 1 "/home/popinet/basilisk-octree/src/grid/cartesian-common.h"
#line 1 "grid/events.h"
#line 1 "/home/popinet/basilisk-octree/src/grid/events.h"





static void event_error (Event * ev, const char * s)
{
  fprintf (qstderr(), "%s:%d: error: %s\n", ev->file, ev->line, s);
  exit (1);
}

static void init_event (Event * ev)
{
  if (ev->arrayi || ev->arrayt) {
    ev->i = ev->t = -1;
    if (ev->arrayi)
      ev->i = ev->arrayi[0];
    else
      ev->t = ev->arrayt[0];
    ev->a = 1;
    ev->expr[1] = NULL;
  }
  else {
    if (ev->nexpr > 0) {
      Expr init = NULL, cond = NULL, inc = NULL;
      for (int j = 0; j < ev->nexpr; j++) {
 int i = -123456; double t = i;
 (* ev->expr[j]) (&i, &t, ev);
 if (i == -123456 && t == -123456) {

   if (cond)
     event_error (ev, "events can only use a single condition");
   cond = ev->expr[j];
 }
 else {

   int i1 = i; double t1 = t;
   (* ev->expr[j]) (&i1, &t1, ev);
   if (i1 == i && t1 == t) {


     if (init)
       event_error (ev, "events can only use a single initialisation");
     init = ev->expr[j];
   }
   else {

     if (inc)
       event_error (ev, "events can only use a single increment");
     inc = ev->expr[j];
   }
 }
      }
      ev->expr[0] = init;
      ev->expr[1] = cond;
      ev->expr[2] = inc;
      ev->nexpr = 0;
    }
    ev->i = ev->t = -1;
    if (ev->expr[0]) {
      (* ev->expr[0]) (&ev->i, &ev->t, ev);
      if (ev->i == 1234567890 || ev->t == 1234567890) {
 ev->i = 1234567890; ev->t = -1;
      }
    }
    else if (ev->expr[2]) {
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (ev->i != -1)
 ev->i = 0;
      if (ev->t != -1)
 ev->t = 0;
    }
  }
}

enum { event_done, event_alive, event_stop };

static int event_finished (Event * ev)
{
  ev->t = ev->i = -1;
  return event_done;
}

void event_register (Event event) {
  assert (Events);
  assert (!event.last);
  int n = 0, found = -1;
  for (Event * ev = Events; !ev->last; ev++) {
    if (found < 0 && !strcmp (event.name, ev->name))
      found = n;
    n++;
  }
  Events = prealloc (Events, (n + 2)*sizeof (Event),__func__,__FILE__,__LINE__);
  Events[n + 1].last = true;
  if (found >= 0)
    for (; n > found; n--)
      Events[n] = Events[n-1];
  Events[n] = event;
  init_event (&Events[n]);
}

static int event_cond (Event * ev, int i, double t)
{
  if (!ev->expr[1])
    return true;
  return (* ev->expr[1]) (&i, &t, ev);
}
#line 120 "/home/popinet/basilisk-octree/src/grid/events.h"
static int event_do (Event * ev, int i, double t, bool action)
{
  if ((i > ev->i && t > ev->t) || !event_cond (ev, i, t))
    return event_finished (ev);
  if (i == ev->i || fabs (t - ev->t) <= 1e-9) {




    if (action && (* ev->action) (i, t, ev)) {
      event_finished (ev);
      return event_stop;
    }
    if (ev->arrayi) {
      ev->i = ev->arrayi[ev->a++];
      if (ev->i < 0)
 return event_finished (ev);
    }
    if (ev->arrayt) {
      ev->t = ev->arrayt[ev->a++];
      if (ev->t < 0)
 return event_finished (ev);
    }
    else if (ev->expr[2]) {
      int i0 = ev->i;
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (i0 == -1 && ev->i != i0)
 ev->i += i + 1;
      if (!event_cond (ev, i + 1, ev->t))
 return event_finished (ev);
    }
    return event_alive;
  }
  return event_alive;
}

static void end_event_do (int i, double t, bool action)
{




  for (Event * ev = Events; !ev->last; ev++)
    if (ev->i == 1234567890 && action) {



      ev->action (i, t, ev);
    }
}

int events (int i, double t, bool action)
{





  if (i == 0)
    for (Event * ev = Events; !ev->last; ev++)
      init_event (ev);

  int inext = 0, cond = 0, cond1 = 0;
  tnext = HUGE;
  for (Event * ev = Events; !ev->last && !cond; ev++)
    if (ev->i != 1234567890 &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond = 1;
  for (Event * ev = Events; !ev->last; ev++) {
    int status = event_do (ev, i, t, action);
    if (status == event_stop) {
      end_event_do (i, t, action);
      return 0;
    }
    if (status == event_alive && ev->i != 1234567890 &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond1 = 1;
    if (ev->t > t && ev->t < tnext)
      tnext = ev->t;
    if (ev->i > i)
      inext = 1;
  }
  if ((!cond || cond1) && (tnext != HUGE || inext))
    return 1;
  end_event_do (i, t, action);
  return 0;
}

void event (const char * name)
{
  for (Event * ev = Events; !ev->last; ev++)
    if (!strcmp (ev->name, name)) {







      (* ev->action) (0, 0, ev);
    }
}

double dtnext (double t, double dt)
{
  if (tnext != HUGE && tnext > t) {
    unsigned int n = (tnext - t)/dt;
    assert (n < INT_MAX);
    if (n == 0)
      dt = tnext - t;
    else {
      double dt1 = (tnext - t)/n;
      if (dt1 > dt + 1e-9)
 dt = (tnext - t)/(n + 1);
      else if (dt1 < dt)
 dt = dt1;
      tnext = t + dt;
    }
  }
  else
    tnext = t + dt;
  return dt;
}
#line 2 "/home/popinet/basilisk-octree/src/grid/cartesian-common.h"

void (* debug) (Point);

#define _val_constant(a,k,l,m) ((const double) _constant[a.i -_NVARMAX])

#undef VARIABLES
#define VARIABLES\
  double Delta = L0*(1./(1 << point.level));\
  \
    double Delta_x = Delta;\
    double Delta_y = Delta;\
    double Delta_z = Delta;\
\
  double x = (ig/2. + (point.i - 2) + 0.5)*Delta + X0; NOT_UNUSED(x);\
\
  double y = (jg/2. + (point.j - 2) + 0.5)*Delta + Y0;\
\
\
\
 NOT_UNUSED(y);\
\
  double z = (kg/2. + (point.k - 2) + 0.5)*Delta + Z0;\
\
\
\
  NOT_UNUSED(z);\
\
  NOT_UNUSED(Delta);\
  \
    NOT_UNUSED(Delta_x);\
    NOT_UNUSED(Delta_y);\
    NOT_UNUSED(Delta_z);\
\
  ;\

#line 32


#line 1 "grid/fpe.h"
#line 1 "/home/popinet/basilisk-octree/src/grid/fpe.h"


#include <signal.h>
#include <unistd.h>

static void gdb()
{
  char command[80];
  sprintf (command, "exec xterm -e gdb -p %d", getpid());
  system (command);
}

static void caught_abort (int sig)
{
  fprintf (qstderr(), "Caught signal %d (Aborted)\n", sig);
  if (last_point.level >= 0) {
    debug (last_point);
    fputc ('\n', qstderr());
  }
  gdb();
}

static void caught_fpe (int sig)
{
  fprintf (qstderr(), "Caught signal %d (Floating Point Exception)\n", sig);
  gdb();
  abort();
}

static void caught_segfault (int sig)
{
  fprintf (qstderr(), "Caught signal %d (Segmentation Fault)\n", sig);
  gdb();
  abort();
}

void catch_fpe (void)
{
  struct sigaction act;
  act.sa_handler = caught_fpe;
  sigemptyset (&act.sa_mask);
  act.sa_flags = 0;
  last_point.level = -1;
  sigaction (8, &act, NULL);
  act.sa_handler = caught_segfault;
  sigaction (11, &act, NULL);
  act.sa_handler = caught_abort;
  act.sa_flags = SA_RESETHAND;
  sigaction (6, &act, NULL);
}
#line 35 "/home/popinet/basilisk-octree/src/grid/cartesian-common.h"

#define end_foreach_face()

scalar new_scalar (const char * name)
{
  int nvar = datasize/sizeof(double);
  scalar s;
  for (s.i = 0; s.i < nvar; s.i++)
    if (!list_lookup (all, s)) {
      all = list_append (all, s);
      init_scalar (s, name);
      trash (((scalar []){s, {-1}}));
      return s;
    }


  assert (nvar < _NVARMAX);
  datasize += sizeof(double); nvar++;
  _attribute = prealloc (_attribute, nvar*sizeof (_Attributes),__func__,__FILE__,__LINE__);
  memset (&_attribute[nvar-1], 0, sizeof (_Attributes));
  all = prealloc (all, sizeof (scalar)*(nvar + 1),__func__,__FILE__,__LINE__);
  s = (scalar){nvar - 1};
  all[nvar - 1] = s;
  all[nvar].i = -1;
  realloc_scalar();
  init_scalar (s, name);
  trash (((scalar []){s, {-1}}));
  return s;
}

scalar new_vertex_scalar (const char * name)
{
  scalar s = new_scalar (name);
  {
#line 68

    _attribute[s.i].d.x = -1;
#line 68

    _attribute[s.i].d.y = -1;
#line 68

    _attribute[s.i].d.z = -1;}
  return s;
}

static vector alloc_vector (const char * name)
{
  vector v;
  char cname[strlen(name) + 3];
  struct { char * x, * y, * z; } ext = {"%s.x", "%s.y", "%s.z"};
  {
#line 78
 {
    sprintf (cname, ext.x, name);
    v.x = new_scalar (cname);
  }
#line 78
 {
    sprintf (cname, ext.y, name);
    v.y = new_scalar (cname);
  }
#line 78
 {
    sprintf (cname, ext.z, name);
    v.z = new_scalar (cname);
  }}
  return v;
}

vector new_vector (const char * name)
{
  vector v = alloc_vector (name);
  init_vector (v, NULL);
  return v;
}

vector new_face_vector (const char * name)
{
  vector v = alloc_vector (name);
  init_face_vector (v, NULL);
  return v;
}

tensor new_tensor (const char * name)
{
  char cname[strlen(name) + 3];
  struct { char * x, * y, * z; } ext = {"%s.x", "%s.y", "%s.z"};
  tensor t;
  {
#line 104
 {
    sprintf (cname, ext.x, name);
    t.x = new_vector (cname);
  }
#line 104
 {
    sprintf (cname, ext.y, name);
    t.y = new_vector (cname);
  }
#line 104
 {
    sprintf (cname, ext.z, name);
    t.z = new_vector (cname);
  }}
  init_tensor (t, NULL);
  return t;
}

tensor new_symmetric_tensor (const char * name)
{
  char cname[strlen(name) + 5];
  struct { char * x, * y, * z; } ext = {"%s.x.x", "%s.y.y", "%s.z.z"};
  tensor t;
  {
#line 117
 {
    sprintf (cname, ext.x, name);
    t.x.x = new_scalar(cname);
  }
#line 117
 {
    sprintf (cname, ext.y, name);
    t.y.y = new_scalar(cname);
  }
#line 117
 {
    sprintf (cname, ext.z, name);
    t.z.z = new_scalar(cname);
  }}

    sprintf (cname, "%s.x.y", name);
    t.x.y = new_scalar(cname);
    t.y.x = t.x.y;


    sprintf (cname, "%s.x.z", name);
    t.x.z = new_scalar(cname);
    t.z.x = t.x.z;
    sprintf (cname, "%s.y.z", name);
    t.y.z = new_scalar(cname);
    t.z.y = t.y.z;

  init_tensor (t, NULL);
  return t;
}

static int nconst = 0;

void init_const_scalar (scalar s, const char * name, double val)
{
  if (s.i - _NVARMAX >= nconst) {
    nconst = s.i - _NVARMAX + 1;
    _constant = prealloc (_constant, nconst*sizeof (double),__func__,__FILE__,__LINE__);
  }
  _constant[s.i - _NVARMAX] = val;
}

scalar new_const_scalar (const char * name, int i, double val)
{
  scalar s = (scalar){i + _NVARMAX};
  init_const_scalar (s, name, val);
  return s;
}

void init_const_vector (vector v, const char * name, double * val)
{
  {
#line 158

    init_const_scalar (v.x, name, *val++);
#line 158

    init_const_scalar (v.y, name, *val++);
#line 158

    init_const_scalar (v.z, name, *val++);}
}

vector new_const_vector (const char * name, int i, double * val)
{
  vector v;
  {
#line 165

    v.x.i = _NVARMAX + i++;
#line 165

    v.y.i = _NVARMAX + i++;
#line 165

    v.z.i = _NVARMAX + i++;}
  init_const_vector (v, name, val);
  return v;
}

void scalar_clone (scalar a, scalar b)
{
  char * name = _attribute[a.i].name;
  double (** boundary) (Point, Point, scalar) = _attribute[a.i].boundary;
  double (** boundary_homogeneous) (Point, Point, scalar) =
    _attribute[a.i].boundary_homogeneous;
  _attribute[a.i] = _attribute[b.i];
  _attribute[a.i].name = name;
  _attribute[a.i].boundary = boundary;
  _attribute[a.i].boundary_homogeneous = boundary_homogeneous;
  for (int i = 0; i < nboundary; i++) {
    _attribute[a.i].boundary[i] = _attribute[b.i].boundary[i];
    _attribute[a.i].boundary_homogeneous[i] = _attribute[b.i].boundary_homogeneous[i];
  }
}

scalar * list_clone (scalar * l)
{
  scalar * list = NULL;
  int nvar = datasize/sizeof(double), map[nvar];
  for (int i = 0; i < nvar; i++)
    map[i] = -1;
  if (l) for (scalar s = *l, *_i27 = l; ((scalar *)&s)->i >= 0; s = *++_i27) {
    scalar c = new_scalar("c");
    scalar_clone (c, s);
    map[s.i] = c.i;
    list = list_append (list, c);
  }
  if (list) for (scalar s = *list, *_i28 = list; ((scalar *)&s)->i >= 0; s = *++_i28)
    {
#line 200

      if (_attribute[s.i].v.x.i >= 0 && map[_attribute[s.i].v.x.i] >= 0)
 _attribute[s.i].v.x.i = map[_attribute[s.i].v.x.i];
#line 200

      if (_attribute[s.i].v.y.i >= 0 && map[_attribute[s.i].v.y.i] >= 0)
 _attribute[s.i].v.y.i = map[_attribute[s.i].v.y.i];
#line 200

      if (_attribute[s.i].v.z.i >= 0 && map[_attribute[s.i].v.z.i] >= 0)
 _attribute[s.i].v.z.i = map[_attribute[s.i].v.z.i];}
  return list;
}

void delete (scalar * list)
{
  if (all == NULL)
    return;

  if (list) for (scalar f = *list, *_i29 = list; ((scalar *)&f)->i >= 0; f = *++_i29) {
    if (_attribute[f.i].delete)
      _attribute[f.i].delete (f);
    pfree (_attribute[f.i].name,__func__,__FILE__,__LINE__); _attribute[f.i].name = NULL;
    pfree (_attribute[f.i].boundary,__func__,__FILE__,__LINE__); _attribute[f.i].boundary = NULL;
    pfree (_attribute[f.i].boundary_homogeneous,__func__,__FILE__,__LINE__); _attribute[f.i].boundary_homogeneous = NULL;
  }

  if (list == all) {
    all[0].i = -1;
    return;
  }

  trash (list);
  if (list) for (scalar f = *list, *_i30 = list; ((scalar *)&f)->i >= 0; f = *++_i30) {
    scalar * s = all;
    for (; s->i >= 0 && s->i != f.i; s++);
    if (s->i == f.i)
      for (; s->i >= 0; s++)
 s[0] = s[1];
  }
}

void free_solver()
{
  delete (all);
  pfree (all,__func__,__FILE__,__LINE__); all = NULL;
  pfree (Events,__func__,__FILE__,__LINE__); Events = NULL;
  pfree (_attribute,__func__,__FILE__,__LINE__); _attribute = NULL;
  pfree (_constant,__func__,__FILE__,__LINE__); _constant = NULL;
  free_grid();
  qpclose_all();
#if TRACE
  trace_off();
#endif
#if MTRACE
  pmuntrace();
#endif
}



void (* boundary_level) (scalar *, int l);
void (* boundary_flux) (vector *);


void boundary (scalar * list)
{ trace ("boundary", "/home/popinet/basilisk-octree/src/grid/cartesian-common.h", 258);
  vector * listf = NULL;
  if (list) for (scalar s = *list, *_i31 = list; ((scalar *)&s)->i >= 0; s = *++_i31)
    if (!is_constant(s) && _attribute[s.i].face)
      listf = vectors_add (listf, _attribute[s.i].v);
  if (listf) {
    boundary_flux (listf);
    pfree (listf,__func__,__FILE__,__LINE__);
  }
  boundary_level (list, -1);
 end_trace("boundary", "/home/popinet/basilisk-octree/src/grid/cartesian-common.h", 268); }

void cartesian_boundary_level (scalar * list, int l)
{
  { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->level) _b->level (_b, list, l); };
}

void cartesian_boundary_flux (vector * list)
{

}

static double symmetry (Point point, Point neighbor, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 281 "/home/popinet/basilisk-octree/src/grid/cartesian-common.h"

  return val(s,0,0,0);
}

static double antisymmetry (Point point, Point neighbor, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 286 "/home/popinet/basilisk-octree/src/grid/cartesian-common.h"

  return -val(s,0,0,0);
}

double (* default_scalar_bc[]) (Point, Point, scalar) = {
  symmetry, symmetry, symmetry, symmetry, symmetry, symmetry
};

scalar cartesian_init_scalar (scalar s, const char * name)
{

  char * pname;
  if (name) {
    pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
    pname = pstrdup (name,__func__,__FILE__,__LINE__);
  }
  else
    pname = _attribute[s.i].name;
  pfree (_attribute[s.i].boundary,__func__,__FILE__,__LINE__);
  pfree (_attribute[s.i].boundary_homogeneous,__func__,__FILE__,__LINE__);

  _attribute[s.i] = (const _Attributes){0};
  _attribute[s.i].name = pname;

  _attribute[s.i].boundary = pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  _attribute[s.i].boundary_homogeneous = pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  for (int b = 0; b < nboundary; b++)
    _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] =
      b < 2*3 ? default_scalar_bc[b] : symmetry;
  _attribute[s.i].gradient = NULL;
  {
#line 316
 {
    _attribute[s.i].d.x = 0;
    _attribute[s.i].v.x.i = -1;
  }
#line 316
 {
    _attribute[s.i].d.y = 0;
    _attribute[s.i].v.y.i = -1;
  }
#line 316
 {
    _attribute[s.i].d.z = 0;
    _attribute[s.i].v.z.i = -1;
  }}
  _attribute[s.i].face = false;
  return s;
}

double (* default_vector_bc[]) (Point, Point, scalar) = {
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry
};

vector cartesian_init_vector (vector v, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  {
#line 333
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      init_scalar (v.x, cname);
    }
    else
      init_scalar (v.x, NULL);
    _attribute[v.x.i].v = v;
  }
#line 333
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      init_scalar (v.y, cname);
    }
    else
      init_scalar (v.y, NULL);
    _attribute[v.y.i].v = v;
  }
#line 333
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.z);
      init_scalar (v.z, cname);
    }
    else
      init_scalar (v.z, NULL);
    _attribute[v.z.i].v = v;
  }}

  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] =
      d < 2*3 ? default_vector_bc[d] : antisymmetry;
  return v;
}

vector cartesian_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
  {
#line 353
 {
    _attribute[v.x.i].d.x = -1;
    _attribute[v.x.i].face = true;
  }
#line 353
 {
    _attribute[v.y.i].d.y = -1;
    _attribute[v.y.i].face = true;
  }
#line 353
 {
    _attribute[v.z.i].d.z = -1;
    _attribute[v.z.i].face = true;
  }}
  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
  return v;
}

tensor cartesian_init_tensor (tensor t, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  {
#line 365
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      init_vector (t.x, cname);
    }
    else
      init_vector (t.x, NULL);
  }
#line 365
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      init_vector (t.y, cname);
    }
    else
      init_vector (t.y, NULL);
  }
#line 365
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.z);
      init_vector (t.z, cname);
    }
    else
      init_vector (t.z, NULL);
  }}
#line 389 "/home/popinet/basilisk-octree/src/grid/cartesian-common.h"
    assert (false);

  return t;
}

void output_cells (FILE * fp)
{
   { foreach(){

#line 396 "/home/popinet/basilisk-octree/src/grid/cartesian-common.h"
 {
    Delta /= 2.;
#line 408 "/home/popinet/basilisk-octree/src/grid/cartesian-common.h"
      for (int i = -1; i <= 1; i += 2) {
 fprintf (fp, "%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n\n",
   x - Delta, y - Delta, z + i*Delta,
   x - Delta, y + Delta, z + i*Delta,
   x + Delta, y + Delta, z + i*Delta,
   x + Delta, y - Delta, z + i*Delta,
   x - Delta, y - Delta, z + i*Delta);
 for (int j = -1; j <= 1; j += 2)
   fprintf (fp, "%g %g %g\n%g %g %g\n\n",
     x + i*Delta, y + j*Delta, z - Delta,
     x + i*Delta, y + j*Delta, z + Delta);
      }

  } } end_foreach(); }
  fflush (fp);
}

void cartesian_debug (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 426 "/home/popinet/basilisk-octree/src/grid/cartesian-common.h"

  char name[80] = "cells";
  if (pid() > 0)
    sprintf (name, "cells-%d", pid());
  FILE * fp = fopen (name, "w");
  output_cells (fp);
  fclose (fp);

  char stencil[80] = "stencil";
  if (pid() > 0)
    sprintf (stencil, "stencil-%d", pid());
  fp = fopen (stencil, "w");
  if (all) for (scalar v = *all, *_i32 = all; ((scalar *)&v)->i >= 0; v = *++_i32)





    fprintf (fp, "x y z %s ", _attribute[v.i].name);

  fputc ('\n', fp);
#line 473 "/home/popinet/basilisk-octree/src/grid/cartesian-common.h"
    for (int k = -2; k <= 2; k++)
      for (int l = -2; l <= 2; l++)
 for (int m = -2; m <= 2; m++) {
   if (all) for (scalar v = *all, *_i33 = all; ((scalar *)&v)->i >= 0; v = *++_i33) {
     fprintf (fp, "%g %g %g ",
       x + k*Delta + _attribute[v.i].d.x*Delta/2.,
       y + l*Delta + _attribute[v.i].d.y*Delta/2.,
       z + m*Delta + _attribute[v.i].d.z*Delta/2.);
     if (allocated(k,l,m))
       fprintf (fp, "%g ", val(v,k,l,m));
     else
       fputs ("n/a ", fp);
   }
   fputc ('\n', fp);
 }

  fclose (fp);

  fp = fopen ("debug.plot", "w");
  fprintf (fp,
    "set term x11\n"
    "set size ratio -1\n"
    "set key outside\n");
  if (all) for (scalar s = *all, *_i34 = all; ((scalar *)&s)->i >= 0; s = *++_i34) {
    char * name = pstrdup (_attribute[s.i].name,__func__,__FILE__,__LINE__), * c = name;
    while (*c != '\0') {
      if (*c == '.')
 *c = '_';
      c++;
    }
    fprintf (fp, "%s = %d\n", name, s.i);
    pfree (name,__func__,__FILE__,__LINE__);
  }
  fclose (fp);

  fprintf (qstderr(),
    "Last point stencils can be displayed using (in gnuplot)\n"
    "  load 'debug.plot'\n"
    "  v=%s\n"







    "  splot '%s' w l lc 0, "
    "'%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 1"
           " title columnhead(4+4*v)",

    _attribute[0].name, name, stencil);
  fflush (qstderr());
}

void cartesian_methods()
{
  init_scalar = cartesian_init_scalar;
  init_vector = cartesian_init_vector;
  init_tensor = cartesian_init_tensor;
  init_face_vector = cartesian_init_face_vector;
  boundary_level = cartesian_boundary_level;
  boundary_flux = cartesian_boundary_flux;
  debug = cartesian_debug;
}

struct _interpolate {
  scalar v;
  double x, y, z;
};


double interpolate (struct _interpolate p)
{ trace ("interpolate", "/home/popinet/basilisk-octree/src/grid/cartesian-common.h", 545);
  Point point = locate ((struct _locate){p.x, p.y, p.z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 546 "/home/popinet/basilisk-octree/src/grid/cartesian-common.h"

  if (point.level < 0)
    { double _ret =  nodata; end_trace("interpolate", "/home/popinet/basilisk-octree/src/grid/cartesian-common.h", 548);  return _ret; }
  scalar v = p.v;
#line 565 "/home/popinet/basilisk-octree/src/grid/cartesian-common.h"
    x = (p.x - x)/Delta - _attribute[v.i].d.x/2.;
    y = (p.y - y)/Delta - _attribute[v.i].d.y/2.;
    z = (p.z - z)/Delta - _attribute[v.i].d.z/2.;
    int i = sign(x), j = sign(y), k = sign(z);
    x = fabs(x); y = fabs(y); z = fabs(z);

    { double _ret =  (((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +
      (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y)*(1. - z) +
     ((val(v,0,0,k)*(1. - x) + val(v,i,0,k)*x)*(1. - y) +
      (val(v,0,j,k)*(1. - x) + val(v,i,j,k)*x)*y)*z); end_trace("interpolate", "/home/popinet/basilisk-octree/src/grid/cartesian-common.h", 574);  return _ret; }

 end_trace("interpolate", "/home/popinet/basilisk-octree/src/grid/cartesian-common.h", 576); }



typedef int bid;

bid new_bid()
{
  int b = nboundary++;
  if (all) for (scalar s = *all, *_i35 = all; ((scalar *)&s)->i >= 0; s = *++_i35) {
    _attribute[s.i].boundary = prealloc (_attribute[s.i].boundary, nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
    _attribute[s.i].boundary_homogeneous = prealloc (_attribute[s.i].boundary_homogeneous,
          nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  }
  if (all) for (scalar s = *all, *_i36 = all; ((scalar *)&s)->i >= 0; s = *++_i36) {
    if (_attribute[s.i].v.x.i < 0)
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] = symmetry;
    else if (_attribute[s.i].v.x.i == s.i) {
      vector v = _attribute[s.i].v;
      {
#line 595

 _attribute[v.y.i].boundary[b] = _attribute[v.y.i].boundary_homogeneous[b] = symmetry;
#line 595

 _attribute[v.z.i].boundary[b] = _attribute[v.z.i].boundary_homogeneous[b] = symmetry;
#line 595

 _attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] = symmetry;}
      _attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] =
 _attribute[v.x.i].face ? NULL : antisymmetry;
    }
  }
  return b;
}



static double periodic_bc (Point point, Point neighbor, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 607 "/home/popinet/basilisk-octree/src/grid/cartesian-common.h"

  return HUGE;
}

void periodic (int dir)
{





    assert (dir <= back);


  int c = dir/2;

  if (all) for (scalar s = *all, *_i37 = all; ((scalar *)&s)->i >= 0; s = *++_i37)
    _attribute[s.i].boundary[2*c] = _attribute[s.i].boundary[2*c + 1] =
      _attribute[s.i].boundary_homogeneous[2*c] = _attribute[s.i].boundary_homogeneous[2*c + 1] =
      periodic_bc;

  if (all) for (scalar s = *all, *_i38 = all; ((scalar *)&s)->i >= 0; s = *++_i38)
    if (_attribute[s.i].face) {
      vector v = _attribute[s.i].v;
      _attribute[v.x.i].boundary[2*c] = _attribute[v.x.i].boundary[2*c + 1] =
 _attribute[v.x.i].boundary_homogeneous[2*c] = _attribute[v.x.i].boundary_homogeneous[2*c + 1] = NULL;
    }

  default_scalar_bc[2*c] = default_scalar_bc[2*c + 1] = periodic_bc;
  default_vector_bc[2*c] = default_vector_bc[2*c + 1] = periodic_bc;
}
#line 4 "/home/popinet/basilisk-octree/src/grid/multigrid-common.h"

#ifndef foreach_level_or_leaf
# define foreach_level_or_leaf foreach_level
# define end_foreach_level_or_leaf end_foreach_level
#endif

#ifndef foreach_coarse_level
# define foreach_coarse_level foreach_level
# define end_foreach_coarse_level end_foreach_level
#endif










static inline void coarsen_average (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 25 "/home/popinet/basilisk-octree/src/grid/multigrid-common.h"

  double sum = 0.;
   { foreach_child()
    sum += val(s,0,0,0); end_foreach_child(); }
  val(s,0,0,0) = sum/(1 << 3);
}

static inline void coarsen_volume_average (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 33 "/home/popinet/basilisk-octree/src/grid/multigrid-common.h"

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 33

  double sum = 0.;
   { foreach_child()
    sum += val_cm(cm,0,0,0)*val(s,0,0,0); end_foreach_child(); }
  val(s,0,0,0) = sum/(1 << 3)/val_cm(cm,0,0,0);
 }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 33

  double sum = 0.;
   { foreach_child()
    sum += val_cm(cm,0,0,0)*val(s,0,0,0); end_foreach_child(); }
  val(s,0,0,0) = sum/(1 << 3)/val_cm(cm,0,0,0);
 }}

static inline void face_average (Point point, vector v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 41 "/home/popinet/basilisk-octree/src/grid/multigrid-common.h"

  {
#line 42
 {







      val(v.x,0,0,0) = (fine(v.x,0,0,0) + fine(v.x,0,1,0) +
        fine(v.x,0,0,1) + fine(v.x,0,1,1))/4.;
      val(v.x,1,0,0) = (fine(v.x,2,0,0) + fine(v.x,2,1,0) +
  fine(v.x,2,0,1) + fine(v.x,2,1,1))/4.;

  }
#line 42
 {







      val(v.y,0,0,0) = (fine(v.y,0,0,0) + fine(v.y,0,0,1) +
        fine(v.y,1,0,0) + fine(v.y,1,0,1))/4.;
      val(v.y,0,1,0) = (fine(v.y,0,2,0) + fine(v.y,0,2,1) +
  fine(v.y,1,2,0) + fine(v.y,1,2,1))/4.;

  }
#line 42
 {







      val(v.z,0,0,0) = (fine(v.z,0,0,0) + fine(v.z,1,0,0) +
        fine(v.z,0,1,0) + fine(v.z,1,1,0))/4.;
      val(v.z,0,0,1) = (fine(v.z,0,0,2) + fine(v.z,1,0,2) +
  fine(v.z,0,1,2) + fine(v.z,1,1,2))/4.;

  }}
}

static inline void coarsen_face (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 59 "/home/popinet/basilisk-octree/src/grid/multigrid-common.h"

  face_average (point, _attribute[s.i].v);
}

static inline void no_coarsen (Point point, scalar s) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 62 "/home/popinet/basilisk-octree/src/grid/multigrid-common.h"
}

static inline void no_data (Point point, scalar s) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 65 "/home/popinet/basilisk-octree/src/grid/multigrid-common.h"

   { foreach_child()
    val(s,0,0,0) = nodata; end_foreach_child(); }
}

void restriction (scalar * list)
{
  scalar * listc = NULL;
  vector * listf = NULL;
  if (list) for (scalar s = *list, *_i39 = list; ((scalar *)&s)->i >= 0; s = *++_i39)
    if (!is_constant(s)) {
      if (_attribute[s.i].face)
 listf = vectors_add (listf, _attribute[s.i].v);
      else
 listc = list_add (listc, s);
    }
  if (listf)
    boundary_flux (listf);
  if (listf || listc) {
    { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->restriction) _b->restriction (_b, list, depth()); };
    for (int l = depth() - 1; l >= 0; l--) {
       { foreach_coarse_level(l){

#line 86 "/home/popinet/basilisk-octree/src/grid/multigrid-common.h"
 {

 if (listc) for (scalar s = *listc, *_i40 = listc; ((scalar *)&s)->i >= 0; s = *++_i40)
   coarsen_average (point, s);
 if (listf) for (vector v = *listf, *_i41 = listf; ((scalar *)&v)->i >= 0; v = *++_i41)
   face_average (point, v);
      } } end_foreach_coarse_level(); }
      { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->restriction) _b->restriction (_b, list, l); };
    }
  }
  pfree (listc,__func__,__FILE__,__LINE__);
  pfree (listf,__func__,__FILE__,__LINE__);
}

void wavelet (scalar s, scalar w)
{
  restriction (((scalar []){s,{-1}}));
   { foreach_fine_to_coarse(){

#line 103 "/home/popinet/basilisk-octree/src/grid/multigrid-common.h"
 {
    double sc[1 << 3];
    int c = 0;
     { foreach_child()
      sc[c++] = val(s,0,0,0); end_foreach_child(); }
    _attribute[s.i].prolongation (point, s);
    c = 0;
     { foreach_child() {

      val(w,0,0,0) = sc[c] - val(s,0,0,0);
      val(s,0,0,0) = sc[c++];
    } end_foreach_child(); }
  } } end_foreach_fine_to_coarse(); }

   { foreach_level(0){

#line 117 "/home/popinet/basilisk-octree/src/grid/multigrid-common.h"
 val(w,0,0,0) = 0.; } end_foreach_level(); }
}

static inline double bilinear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 121 "/home/popinet/basilisk-octree/src/grid/multigrid-common.h"








    return (27.*coarse(s,0,0,0) +
     9.*(coarse(s,child.x,0,0) + coarse(s,0,child.y,0) +
  coarse(s,0,0,child.z)) +
     3.*(coarse(s,child.x,child.y,0) + coarse(s,child.x,0,child.z) +
  coarse(s,0,child.y,child.z)) +
     coarse(s,child.x,child.y,child.z))/64.;

}

static inline void refine_bilinear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 139 "/home/popinet/basilisk-octree/src/grid/multigrid-common.h"

   { foreach_child()
    val(s,0,0,0) = bilinear (point, s); end_foreach_child(); }
}

static inline double biquadratic (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 145 "/home/popinet/basilisk-octree/src/grid/multigrid-common.h"

#line 156 "/home/popinet/basilisk-octree/src/grid/multigrid-common.h"
  assert (false);
  return 0.;

}

static inline void refine_biquadratic (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 162 "/home/popinet/basilisk-octree/src/grid/multigrid-common.h"

   { foreach_child()
    val(s,0,0,0) = biquadratic (point, s); end_foreach_child(); }
}

static inline void refine_linear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 168 "/home/popinet/basilisk-octree/src/grid/multigrid-common.h"

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 168

  coord g;
  if (_attribute[s.i].gradient)
    {
#line 171

      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
#line 171

      g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));
#line 171

      g.z = _attribute[s.i].gradient (val(s,0,0,-1), val(s,0,0,0), val(s,0,0,1));}
  else
    {
#line 174

      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
#line 174

      g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;
#line 174

      g.z = (val(s,0,0,1) - val(s,0,0,-1))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*val_cm(cm,0,0,0), sum = val_cm(cm,0,0,0)*(1 << 3);
   { foreach_child() {
    val(s,0,0,0) = sc;
    {
#line 180

      val(s,0,0,0) += child.x*g.x*val_cm(cm,-child.x,0,0)/cmc;
#line 180

      val(s,0,0,0) += child.y*g.y*val_cm(cm,0,-child.y,0)/cmc;
#line 180

      val(s,0,0,0) += child.z*g.z*val_cm(cm,0,0,-child.z)/cmc;}
    sum -= val_cm(cm,0,0,0);
  } end_foreach_child(); }
  assert (fabs(sum) < 1e-10);
 }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 168

  coord g;
  if (_attribute[s.i].gradient)
    {
#line 171

      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
#line 171

      g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));
#line 171

      g.z = _attribute[s.i].gradient (val(s,0,0,-1), val(s,0,0,0), val(s,0,0,1));}
  else
    {
#line 174

      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
#line 174

      g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;
#line 174

      g.z = (val(s,0,0,1) - val(s,0,0,-1))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*val_cm(cm,0,0,0), sum = val_cm(cm,0,0,0)*(1 << 3);
   { foreach_child() {
    val(s,0,0,0) = sc;
    {
#line 180

      val(s,0,0,0) += child.x*g.x*val_cm(cm,-child.x,0,0)/cmc;
#line 180

      val(s,0,0,0) += child.y*g.y*val_cm(cm,0,-child.y,0)/cmc;
#line 180

      val(s,0,0,0) += child.z*g.z*val_cm(cm,0,0,-child.z)/cmc;}
    sum -= val_cm(cm,0,0,0);
  } end_foreach_child(); }
  assert (fabs(sum) < 1e-10);
 }}

static inline void refine_reset (Point point, scalar v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 188 "/home/popinet/basilisk-octree/src/grid/multigrid-common.h"

   { foreach_child()
    val(v,0,0,0) = 0.; end_foreach_child(); }
}

static inline void refine_injection (Point point, scalar v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 194 "/home/popinet/basilisk-octree/src/grid/multigrid-common.h"

  double val = val(v,0,0,0);
   { foreach_child()
    val(v,0,0,0) = val; end_foreach_child(); }
}

vector multigrid_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  {
#line 203

    _attribute[v.y.i].coarsen = no_coarsen;
#line 203

    _attribute[v.z.i].coarsen = no_coarsen;
#line 203

    _attribute[v.x.i].coarsen = no_coarsen;}
  _attribute[v.x.i].coarsen = coarsen_face;
  return v;
}

void multigrid_debug (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 210 "/home/popinet/basilisk-octree/src/grid/multigrid-common.h"

  cartesian_debug (point);

  if (point.level > 0) {
    char name[80] = "coarse";
    if (pid() > 0)
      sprintf (name, "coarse-%d", pid());
    FILE * fp = fopen (name, "w");
#line 240 "/home/popinet/basilisk-octree/src/grid/multigrid-common.h"
      double xc = x - child.x*Delta/2., yc = y - child.y*Delta/2.;
      double zc = z - child.z*Delta/2.;
      for (int k = 0; k <= 1; k++)
 for (int l = 0; l <= 1; l++)
   for (int m = 0; m <= 1; m++) {
     if (all) for (scalar v = *all, *_i42 = all; ((scalar *)&v)->i >= 0; v = *++_i42)
       fprintf (fp, "%g %g %g %g ",
         xc + k*child.x*Delta*2. + _attribute[v.i].d.x*Delta,
         yc + l*child.y*Delta*2. + _attribute[v.i].d.y*Delta,
         zc + m*child.z*Delta*2. + _attribute[v.i].d.z*Delta,
         coarse(v,k*child.x,l*child.y,m*child.z));
     fputc ('\n', fp);
   }
      fprintf (qstderr(), ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 3 t ''",
        name);

    fclose (fp);
  }

  if (is_coarse()) {
    char name[80] = "fine";
    if (pid() > 0)
      sprintf (name, "fine-%d", pid());
    FILE * fp = fopen (name, "w");
#line 293 "/home/popinet/basilisk-octree/src/grid/multigrid-common.h"
      double xf = x - Delta/4., yf = y - Delta/4., zf = z - Delta/4.;
      for (int k = -2; k <= 3; k++)
 for (int l = -2; l <= 3; l++)
   for (int m = -2; m <= 3; m++) {
     if (all) for (scalar v = *all, *_i43 = all; ((scalar *)&v)->i >= 0; v = *++_i43) {
       fprintf (fp, "%g %g %g ",
         xf + k*Delta/2. + _attribute[v.i].d.x*Delta/4.,
         yf + l*Delta/2. + _attribute[v.i].d.y*Delta/4.,
         zf + m*Delta/2. + _attribute[v.i].d.z*Delta/4.);
       if (allocated_child(k,l,m))
  fprintf (fp, "%g ", fine(v,k,l,m));
       else
  fputs ("n/a ", fp);
     }
     fputc ('\n', fp);
   }
      fprintf (qstderr(), ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 2 t ''",
        name);

    fclose (fp);
  }
  fflush (qstderr());
}

void multigrid_methods()
{
  cartesian_methods();
  debug = multigrid_debug;
  init_face_vector = multigrid_init_face_vector;
}
#line 4 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"












int refine_cell (Point point, scalar * list, int flag, Cache * refined)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 17 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"

  int nr = 0;


  if (level > 0)
    for (int k = 0; k != 2*child.x; k += child.x)

      for (int l = 0; l != 2*child.y; l += child.y)


 for (int m = 0; m != 2*child.z; m += child.z)

   if (aparent(k,l,m).pid >= 0 && is_leaf(aparent(k,l,m))) {
     Point p = point;


     p.level = point.level - 1;
     p.i = (point.i + 2)/2 + k;

       p.j = (point.j + 2)/2 + l;


       p.k = (point.k + 2)/2 + m;

       nr += refine_cell (p, list, flag, refined);
     aparent(k,l,m).flags |= flag;
   }



  cell.flags &= ~leaf;


  increment_neighbors (point);

  int cflag = is_active(cell) ? (active|leaf) : leaf;
   { foreach_child()
    cell.flags |= cflag; end_foreach_child(); }


  if (list) for (scalar s = *list, *_i44 = list; ((scalar *)&s)->i >= 0; s = *++_i44)
    if (is_local(cell) || _attribute[s.i].face)
      _attribute[s.i].refine (point, s);

#if _MPI
  if (is_border(cell)) {
     { foreach_child() {
      bool bord = false;
       { foreach_neighbor() {
 if (!is_local(cell) || (level > 0 && !is_local(aparent(0,0,0))))
   bord = true, foreach_neighbor_break();
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
    { foreach_child()
     if (!is_local(cell))
       bord = true, foreach_child_break(); end_foreach_child(); }
 if (bord)
   foreach_neighbor_break();
      } end_foreach_neighbor(); }
      if (bord)
 cell.flags |= border;
    } end_foreach_child(); }
    if (refined)
      cache_append (refined, point, cell.flags);
    nr++;
  }
#endif
  return nr;
}

bool coarsen_cell (Point point, scalar * list)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 87 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"




  int pid = cell.pid;
   { foreach_child()
    if (cell.neighbors || (cell.pid < 0 && cell.pid != pid))
      return false; end_foreach_child(); }



  if (list) for (scalar s = *list, *_i45 = list; ((scalar *)&s)->i >= 0; s = *++_i45)
    _attribute[s.i].coarsen (point, s);


  cell.flags |= leaf;
  cell.flags &= ~halo;


  decrement_neighbors (point);

#if _MPI
  if (!is_local(cell)) {
    cell.flags &= ~(active|border);
     { foreach_neighbor(1)
      if (cell.neighbors)
  { foreach_child()
   if (allocated(0,0,0) && is_local(cell) && !is_border(cell))
     cell.flags |= border; end_foreach_child(); } end_foreach_neighbor(); }
  }
#endif

  return true;
}

void coarsen_cell_recursive (Point point, scalar * list)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 123 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"



   { foreach_child()
    if (cell.neighbors)
       { foreach_neighbor(1)
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
   coarsen_cell_recursive (point, list); end_foreach_neighbor(); } end_foreach_child(); }

  assert (coarsen_cell (point, list));
}

#if _MPI
void mpi_boundary_refine (scalar *, int);
void mpi_boundary_coarsen (int, int);
void mpi_boundary_update (void);
bool balance();
#else
# define mpi_boundary_refine(...)
# define mpi_boundary_coarsen(...)
# define mpi_boundary_update()
# define balance(...) false
#endif

typedef struct {
  int nc, nf;
} astats;

struct Adapt {
  scalar * slist;
  double * max;
  int maxlevel;
  int minlevel;
  scalar * list;
  scalar * listb;
};

astats adapt_wavelet (struct Adapt p)
{
  scalar * listcm = NULL;

  if (is_constant(cm))
    restriction (p.slist);
  else {
    if (p.list == NULL)
      listcm = list_concat (NULL, ((scalar []){cm,fm.x,fm.y,fm.z,{-1}}));
    scalar * listr = list_concat (p.slist, ((scalar []){cm,{-1}}));
    restriction (listr);
    pfree (listr,__func__,__FILE__,__LINE__);
  }
  if (p.list == NULL) {
    if (all) for (scalar s = *all, *_i46 = all; ((scalar *)&s)->i >= 0; s = *++_i46)
      listcm = list_add (listcm, s);
    p.list = listcm;
  }

  astats st = {0, 0};
  scalar * listc = NULL;
  if (p.list) for (scalar s = *p.list, *_i47 = p.list; ((scalar *)&s)->i >= 0; s = *++_i47)
    if (!is_constant(s) && _attribute[s.i].coarsen != no_coarsen)
      listc = list_add (listc, s);


  if (p.minlevel < 1)
    p.minlevel = 1;
  ((Quadtree *)grid)->refined.n = 0;
  static const int refined = 1 << user, too_fine = 1 << (user + 1);
   { foreach_cell(){

#line 190 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"
 {
    if (is_active(cell)) {
      static const int too_coarse = 1 << (user + 2);
      if (is_leaf (cell)) {
 if (cell.flags & too_coarse) {
   cell.flags &= ~too_coarse;
   refine_cell (point, listc, refined, &((Quadtree *)grid)->refined);
   st.nf++;
 }
 continue;
      }
      else {
 if (cell.flags & refined) {

   cell.flags &= ~too_coarse;
   continue;
 }

 bool local = is_local(cell);
 if (!local)
    { foreach_child()
     if (is_local(cell))
       local = true, foreach_child_break(); end_foreach_child(); }
 if (local) {
   int i = 0;
   static const int just_fine = 1 << (user + 3);
   if (p.slist) for (scalar s = *p.slist, *_i48 = p.slist; ((scalar *)&s)->i >= 0; s = *++_i48) {
     double max = p.max[i++], sc[1 << 3];
     int c = 0;
      { foreach_child()
       sc[c++] = val(s,0,0,0); end_foreach_child(); }
     _attribute[s.i].prolongation (point, s);
     c = 0;
      { foreach_child() {
       double e = fabs(sc[c] - val(s,0,0,0));
       if (e > max) {
  cell.flags &= ~too_fine;
  cell.flags |= too_coarse;
       }
       else if (e <= max/1.5 && !(cell.flags & (too_coarse|just_fine))) {
  if (level >= p.minlevel)
    cell.flags |= too_fine;
       }
       else if (!(cell.flags & too_coarse)) {
  cell.flags &= ~too_fine;
  cell.flags |= just_fine;
       }
       val(s,0,0,0) = sc[c++];
     } end_foreach_child(); }
   }
    { foreach_child() {
     cell.flags &= ~just_fine;
     if (!is_leaf(cell) || !is_active(cell) || level == p.maxlevel)
       cell.flags &= ~too_coarse;
   } end_foreach_child(); }
 }

 if (level == p.maxlevel - 1)
   continue;
      }
    }
    else
      continue;
  } } end_foreach_cell(); }
  mpi_boundary_refine (listc, 1);



  for (int l = depth(); l >= p.minlevel; l--) {
     { foreach_cell(){

#line 259 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"

      if (!(cell.pid < 0)) {
 if (level == l) {
   if (!is_leaf(cell)) {
     if (cell.flags & refined)

       cell.flags &= ~(refined|too_fine);
     else if (cell.flags & too_fine) {
       if (is_local(cell) && coarsen_cell (point, listc))
  st.nc++;
       cell.flags &= ~too_fine;
     }
   }
   if (cell.flags & too_fine)
     cell.flags &= ~too_fine;
   else if (aparent(0,0,0).flags & too_fine)
     aparent(0,0,0).flags &= ~too_fine;
   continue;
 }
 else if (is_leaf(cell))
   continue;
      } } end_foreach_cell(); }
    mpi_boundary_coarsen (l, too_fine);
  }
  pfree (listc,__func__,__FILE__,__LINE__);

  mpi_all_reduce (st.nf, MPI_INT, MPI_SUM);
  mpi_all_reduce (st.nc, MPI_INT, MPI_SUM);
  if (st.nc || st.nf) {
    mpi_boundary_update();
    boundary (p.listb ? p.listb : p.list);
    while (balance());
  }
  pfree (listcm,__func__,__FILE__,__LINE__);

  return st;
}
#line 318 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"
static void refine_level (int depth)
{
  do { int refined; do { refined = 0; ((Quadtree *)grid)->refined.n = 0;  { foreach_leaf(){

#line 320 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"
 if (level < depth) { refine_cell (point, NULL, 0, &((Quadtree *)grid)->refined); refined++; continue; } } end_foreach_leaf(); } mpi_all_reduce (refined, MPI_INT, MPI_SUM); if (refined) { mpi_boundary_refine (NULL, 0); mpi_boundary_update(); boundary (NULL); while (balance()); } } while (refined); } while(0);
}
#line 343 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"
static void halo_restriction_flux (vector * list)
{
  vector * listv = NULL;
  if (list) for (vector v = *list, *_i49 = list; ((scalar *)&v)->i >= 0; v = *++_i49)
    if (!is_constant(v.x))
      listv = vectors_add (listv, v);

  if (listv) {
    for (int l = depth() - 1; l >= 0; l--)
       { foreach_halo (prolongation, l){

#line 352 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"

 {
#line 353
 {
#line 369 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"
   if ((!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0))
     if (listv) for (vector f = *listv, *_i50 = listv; ((scalar *)&f)->i >= 0; f = *++_i50)
       val(f.x,0,0,0) = (fine(f.x,0,0,0) + fine(f.x,0,1,0) +
         fine(f.x,0,0,1) + fine(f.x,0,1,1))/4.;
   if ((!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0))
     if (listv) for (vector f = *listv, *_i51 = listv; ((scalar *)&f)->i >= 0; f = *++_i51)
       val(f.x,1,0,0) = (fine(f.x,2,0,0) + fine(f.x,2,1,0) +
   fine(f.x,2,0,1) + fine(f.x,2,1,1))/4.;

      }
#line 353
 {
#line 369 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"
   if ((!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0))
     if (listv) for (vector f = *listv, *_i50 = listv; ((scalar *)&f)->i >= 0; f = *++_i50)
       val(f.y,0,0,0) = (fine(f.y,0,0,0) + fine(f.y,0,0,1) +
         fine(f.y,1,0,0) + fine(f.y,1,0,1))/4.;
   if ((!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0))
     if (listv) for (vector f = *listv, *_i51 = listv; ((scalar *)&f)->i >= 0; f = *++_i51)
       val(f.y,0,1,0) = (fine(f.y,0,2,0) + fine(f.y,0,2,1) +
   fine(f.y,1,2,0) + fine(f.y,1,2,1))/4.;

      }
#line 353
 {
#line 369 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"
   if ((!is_leaf (neighbor(0,0,-1)) && neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0))
     if (listv) for (vector f = *listv, *_i50 = listv; ((scalar *)&f)->i >= 0; f = *++_i50)
       val(f.z,0,0,0) = (fine(f.z,0,0,0) + fine(f.z,1,0,0) +
         fine(f.z,0,1,0) + fine(f.z,1,1,0))/4.;
   if ((!is_leaf (neighbor(0,0,1)) && neighbor(0,0,1).neighbors && neighbor(0,0,1).pid >= 0))
     if (listv) for (vector f = *listv, *_i51 = listv; ((scalar *)&f)->i >= 0; f = *++_i51)
       val(f.z,0,0,1) = (fine(f.z,0,0,2) + fine(f.z,1,0,2) +
   fine(f.z,0,1,2) + fine(f.z,1,1,2))/4.;

      }} } end_foreach_halo(); }
    pfree (listv,__func__,__FILE__,__LINE__);
  }
}



static scalar quadtree_init_scalar (scalar s, const char * name)
{
  s = cartesian_init_scalar (s, name);
  _attribute[s.i].refine = _attribute[s.i].prolongation = refine_bilinear;
  _attribute[s.i].coarsen = coarsen_average;
  return s;
}


#line 393

static void refine_face_x (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 395 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"

  assert (3 > 1);
  vector v = _attribute[s.i].v;
  if (!(!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(-1,0,0)))) {
    double g1 = (val(v.x,0,+1,0) - val(v.x,0,-1,0))/8.;
    double g2 = (val(v.x,0,0,+1) - val(v.x,0,0,-1))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.x,0,j,k) = val(v.x,0,0,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (!(!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0) && neighbor(1,0,0).neighbors &&
      (is_local(cell) || is_local(neighbor(1,0,0)))) {
    double g1 = (val(v.x,1,+1,0) - val(v.x,1,-1,0))/8.;
    double g2 = (val(v.x,1,0,+1) - val(v.x,1,0,-1))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.x,2,j,k) = val(v.x,1,0,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (is_local(cell)) {
    double g1 = (val(v.x,0,+1,0) + val(v.x,1,+1,0) - val(v.x,0,-1,0) - val(v.x,1,-1,0))/16.;
    double g2 = (val(v.x,0,0,+1) + val(v.x,1,0,+1) - val(v.x,0,0,-1) - val(v.x,1,0,-1))/16.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.x,1,j,k) = (val(v.x,0,0,0) + val(v.x,1,0,0))/2. + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
}
#line 393

static void refine_face_y (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 395 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"

  assert (3 > 1);
  vector v = _attribute[s.i].v;
  if (!(!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(0,-1,0)))) {
    double g1 = (val(v.y,0,0,+1) - val(v.y,0,0,-1))/8.;
    double g2 = (val(v.y,+1,0,0) - val(v.y,-1,0,0))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.y,k,0,j) = val(v.y,0,0,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (!(!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0) && neighbor(0,1,0).neighbors &&
      (is_local(cell) || is_local(neighbor(0,1,0)))) {
    double g1 = (val(v.y,0,1,+1) - val(v.y,0,1,-1))/8.;
    double g2 = (val(v.y,+1,1,0) - val(v.y,-1,1,0))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.y,k,2,j) = val(v.y,0,1,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (is_local(cell)) {
    double g1 = (val(v.y,0,0,+1) + val(v.y,0,1,+1) - val(v.y,0,0,-1) - val(v.y,0,1,-1))/16.;
    double g2 = (val(v.y,+1,0,0) + val(v.y,+1,1,0) - val(v.y,-1,0,0) - val(v.y,-1,1,0))/16.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.y,k,1,j) = (val(v.y,0,0,0) + val(v.y,0,1,0))/2. + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
}
#line 393

static void refine_face_z (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 395 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"

  assert (3 > 1);
  vector v = _attribute[s.i].v;
  if (!(!is_leaf (neighbor(0,0,-1)) && neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(0,0,-1)))) {
    double g1 = (val(v.z,+1,0,0) - val(v.z,-1,0,0))/8.;
    double g2 = (val(v.z,0,+1,0) - val(v.z,0,-1,0))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.z,j,k,0) = val(v.z,0,0,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (!(!is_leaf (neighbor(0,0,1)) && neighbor(0,0,1).neighbors && neighbor(0,0,1).pid >= 0) && neighbor(0,0,1).neighbors &&
      (is_local(cell) || is_local(neighbor(0,0,1)))) {
    double g1 = (val(v.z,+1,0,1) - val(v.z,-1,0,1))/8.;
    double g2 = (val(v.z,0,+1,1) - val(v.z,0,-1,1))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.z,j,k,2) = val(v.z,0,0,1) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (is_local(cell)) {
    double g1 = (val(v.z,+1,0,0) + val(v.z,+1,0,1) - val(v.z,-1,0,0) - val(v.z,-1,0,1))/16.;
    double g2 = (val(v.z,0,+1,0) + val(v.z,0,+1,1) - val(v.z,0,-1,0) - val(v.z,0,-1,1))/16.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.z,j,k,1) = (val(v.z,0,0,0) + val(v.z,0,0,1))/2. + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
}

void refine_face (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 424 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"

  vector v = _attribute[s.i].v;
  {
#line 426

    _attribute[v.x.i].prolongation (point, v.x);
#line 426

    _attribute[v.y.i].prolongation (point, v.y);
#line 426

    _attribute[v.z.i].prolongation (point, v.z);}
}

void refine_face_solenoidal (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 431 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"

  refine_face (point, s);

  if (is_local(cell)) {

    vector v = _attribute[s.i].v;
    double d[1 << 3], p[1 << 3];
    int i = 0;
     { foreach_child() {
      d[i] = 0.;
      {
#line 441

 d[i] += val(v.x,1,0,0) - val(v.x,0,0,0);
#line 441

 d[i] += val(v.y,0,1,0) - val(v.y,0,0,0);
#line 441

 d[i] += val(v.z,0,0,1) - val(v.z,0,0,0);}
      i++;
    } end_foreach_child(); }
#line 455 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"
    static double m[7][7] = {
      {7./12,5./24,3./8,5./24,3./8,1./4,1./3},
      {5./24,7./12,3./8,5./24,1./4,3./8,1./3},
      {3./8,3./8,3./4,1./4,3./8,3./8,1./2},
      {5./24,5./24,1./4,7./12,3./8,3./8,1./3},
      {3./8,1./4,3./8,3./8,3./4,3./8,1./2},
      {1./4,3./8,3./8,3./8,3./8,3./4,1./2},
      {1./3,1./3,1./2,1./3,1./2,1./2,5./6}
    };
    p[0] = 0.;
    for (int i = 0; i < 7; i++) {
      p[i + 1] = 0.;
      for (int j = 0; j < 7; j++)
 p[i + 1] += m[i][j]*d[j+1];
    }
    for (int k = 0; k <= 1; k++) {
      fine(v.x,1,0,k) += p[4+k] - p[0+k];
      fine(v.x,1,1,k) += p[6+k] - p[2+k];
      fine(v.y,0,1,k) += p[2+k] - p[0+k];
      fine(v.y,1,1,k) += p[6+k] - p[4+k];
    }
    fine(v.z,0,0,1) += p[1] - p[0];
    fine(v.z,0,1,1) += p[3] - p[2];
    fine(v.z,1,0,1) += p[5] - p[4];
    fine(v.z,1,1,1) += p[7] - p[6];

  }

}

vector quadtree_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  {
#line 488

    _attribute[v.x.i].coarsen = _attribute[v.x.i].refine = no_coarsen;
#line 488

    _attribute[v.y.i].coarsen = _attribute[v.y.i].refine = no_coarsen;
#line 488

    _attribute[v.z.i].coarsen = _attribute[v.z.i].refine = no_coarsen;}
  _attribute[v.x.i].coarsen = coarsen_face;
  _attribute[v.x.i].refine = refine_face;
  {
#line 492

    _attribute[v.x.i].prolongation = refine_face_x;
#line 492

    _attribute[v.y.i].prolongation = refine_face_y;
#line 492

    _attribute[v.z.i].prolongation = refine_face_z;}
  return v;
}

static void quadtree_boundary_level (scalar * list, int l)
{
  int depth = l < 0 ? depth() : l;

  scalar * listdef = NULL, * listc = NULL, * list2 = NULL;
  if (list) for (scalar s = *list, *_i52 = list; ((scalar *)&s)->i >= 0; s = *++_i52)
    if (!is_constant (s)) {
      if (_attribute[s.i].coarsen == coarsen_average) {
 listdef = list_add (listdef, s);
 list2 = list_add (list2, s);
      }
      else if (_attribute[s.i].coarsen != no_coarsen) {
 listc = list_add (listc, s);
 if (_attribute[s.i].face)
   {
#line 511

     list2 = list_add (list2, _attribute[s.i].v.x);
#line 511

     list2 = list_add (list2, _attribute[s.i].v.y);
#line 511

     list2 = list_add (list2, _attribute[s.i].v.z);}
 else
   list2 = list_add (list2, s);
      }
    }

  if (listdef || listc) {
    { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->halo_restriction) _b->halo_restriction (_b, list2, depth); };
    for (int i = depth - 1; i >= 0; i--) {
       { foreach_halo (restriction, i){

#line 521 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"
 {
 if (listdef) for (scalar s = *listdef, *_i53 = listdef; ((scalar *)&s)->i >= 0; s = *++_i53)
   coarsen_average (point, s);
 if (listc) for (scalar s = *listc, *_i54 = listc; ((scalar *)&s)->i >= 0; s = *++_i54)
   _attribute[s.i].coarsen (point, s);
      } } end_foreach_halo(); }
      { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->halo_restriction) _b->halo_restriction (_b, list2, i); };
    }
    pfree (listdef,__func__,__FILE__,__LINE__);
    pfree (listc,__func__,__FILE__,__LINE__);
    pfree (list2,__func__,__FILE__,__LINE__);
  }

  scalar * listr = NULL;
  vector * listf = NULL;
  if (list) for (scalar s = *list, *_i55 = list; ((scalar *)&s)->i >= 0; s = *++_i55)
    if (!is_constant (s) && _attribute[s.i].refine != no_coarsen) {
      if (_attribute[s.i].face)
 listf = vectors_add (listf, _attribute[s.i].v);
      else
 listr = list_add (listr, s);
    }

  if (listr || listf) {
    { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->halo_prolongation) _b->halo_prolongation (_b, list, 0, l); };
    for (int i = 0; i < depth; i++) {
       { foreach_halo (prolongation, i){

#line 547 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"
 {
 if (listr) for (scalar s = *listr, *_i56 = listr; ((scalar *)&s)->i >= 0; s = *++_i56)
          _attribute[s.i].prolongation (point, s);
 if (listf) for (vector v = *listf, *_i57 = listf; ((scalar *)&v)->i >= 0; v = *++_i57)
   {
#line 551

     _attribute[v.x.i].prolongation (point, v.x);
#line 551

     _attribute[v.y.i].prolongation (point, v.y);
#line 551

     _attribute[v.z.i].prolongation (point, v.z);}
      } } end_foreach_halo(); }
      { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->halo_prolongation) _b->halo_prolongation (_b, list, i + 1, l); };
    }
    pfree (listr,__func__,__FILE__,__LINE__);
    pfree (listf,__func__,__FILE__,__LINE__);
  }
}

double treex (Point point) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 561 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"

  if (level == 0)
    return 0;
  double i = 2*child.x - child.y;
  if (i <= 1 && i >= -1) i = -i;
  return treex(parent) + i/(1 << 2*(level - 1));
}

double treey (Point point) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 569 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"

  if (level == 0)
    return 0;
  return treey(parent) + 4./(1 << 2*(level - 1));
}

void output_tree (FILE * fp)
{
   { foreach_cell(){

#line 577 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"

    if (cell.neighbors)
       { foreach_child()
 if (is_local(cell))
   fprintf (fp, "%g %g\n%g %g\n\n",
     treex(parent), treey(parent), treex(point), treey(point)); end_foreach_child(); }; } end_foreach_cell(); }
}

void quadtree_check()
{


  long nleaves = 0, nactive = 0;
   { foreach_cell_all(){

#line 590 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"
 {
    if (is_leaf(cell)) {
      assert (cell.pid >= 0);
      nleaves++;
    }
    if (is_local(cell))
      assert (is_active(cell) || (!is_leaf(cell) && !cell.neighbors && cell.pid >= 0));
    if (is_active(cell))
      nactive++;

    int neighbors = 0;
     { foreach_neighbor(1)
      if (allocated(0,0,0) && (!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
 neighbors++; end_foreach_neighbor(); }
    assert (cell.neighbors == neighbors);


    if (!cell.neighbors)
      assert (!allocated_child(0,0,0));
  } } end_foreach_cell_all(); }


  long reachable = 0;
   { foreach_cell(){

#line 613 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"
 {
    if (is_active(cell))
      reachable++;
    else
      continue;
  } } end_foreach_cell(); }
  assert (nactive == reachable);


  reachable = 0;
   { foreach_cell(){

#line 623 "/home/popinet/basilisk-octree/src/grid/quadtree-common.h"

    if (is_leaf(cell)) {
      reachable++;
      continue;
    } } end_foreach_cell(); }
  assert (nleaves == reachable);
}

void quadtree_methods()
{
  multigrid_methods();
  init_scalar = quadtree_init_scalar;
  init_face_vector = quadtree_init_face_vector;
  boundary_level = quadtree_boundary_level;
  boundary_flux = halo_restriction_flux;
}
#line 1713 "/home/popinet/basilisk-octree/src/grid/tree.h"

#if _MPI
#line 1 "grid/quadtree-mpi.h"
#line 1 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"

int debug_iteration = -1;

void debug_mpi (FILE * fp1);

typedef struct {
  CacheLevel * halo;
  void * buf;
  MPI_Request r;
  int depth;
  int pid;
  int maxdepth;
} Rcv;

typedef struct {
  Rcv * rcv;
  char * name;
  int npid;
} RcvPid;

typedef struct {
  RcvPid * rcv, * snd;
} SndRcv;

typedef struct {
  Boundary parent;

  SndRcv restriction, restriction_root, halo_restriction, prolongation;
  Array * send, * receive;
} MpiBoundary;

static void cache_level_init (CacheLevel * c)
{
  c->p = NULL;
  c->n = c->nm = 0;
}

static void rcv_append (Point point, Rcv * rcv)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 39 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"

  if (level > rcv->depth) {
    rcv->halo = prealloc (rcv->halo, (level + 1)*sizeof (CacheLevel),__func__,__FILE__,__LINE__);
    for (int j = rcv->depth + 1; j <= level; j++)
      cache_level_init (&rcv->halo[j]);
    rcv->depth = level;
  }
  cache_level_append (&rcv->halo[level], point);
  if (level > rcv->maxdepth)
    rcv->maxdepth = level;
}

void rcv_print (Rcv * rcv, FILE * fp, const char * prefix)
{
  for (int l = 0; l <= rcv->depth; l++)
    if (rcv->halo[l].n > 0)
       { foreach_cache_level(rcv->halo[l], l,){

#line 55 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"

 fprintf (fp, "%s%g %g %g %d %d\n", prefix, x, y, z, rcv->pid, level); } end_foreach_cache_level(); }
}

static void rcv_free_buf (Rcv * rcv)
{
  if (rcv->buf) {
    prof_start ("rcv_pid_receive");
    MPI_Wait (&rcv->r, MPI_STATUS_IGNORE);
    pfree (rcv->buf,__func__,__FILE__,__LINE__);
    rcv->buf = NULL;
    prof_stop();
  }
}

static void rcv_destroy (Rcv * rcv)
{
  rcv_free_buf (rcv);
  for (int i = 0; i <= rcv->depth; i++)
    if (rcv->halo[i].n > 0)
      pfree (rcv->halo[i].p,__func__,__FILE__,__LINE__);
  pfree (rcv->halo,__func__,__FILE__,__LINE__);
}

static RcvPid * rcv_pid_new (const char * name)
{
  RcvPid * r = pcalloc (sizeof (RcvPid), 1,__func__,__FILE__,__LINE__);
  r->name = pstrdup (name,__func__,__FILE__,__LINE__);
  return r;
}

static Rcv * rcv_pid_pointer (RcvPid * p, int pid)
{
  assert (pid >= 0 && pid < npe());

  int i;
  for (i = 0; i < p->npid; i++)
    if (pid == p->rcv[i].pid)
      break;

  if (i == p->npid) {
    p->rcv = prealloc (p->rcv, ++p->npid*sizeof (Rcv),__func__,__FILE__,__LINE__);
    Rcv * rcv = &p->rcv[p->npid-1];
    rcv->pid = pid;
    rcv->depth = rcv->maxdepth = 0;
    rcv->halo = pmalloc (sizeof (CacheLevel),__func__,__FILE__,__LINE__);
    rcv->buf = NULL;
    cache_level_init (&rcv->halo[0]);
  }
  return &p->rcv[i];
}

static void rcv_pid_append (RcvPid * p, int pid, Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 108 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"

  rcv_append (point, rcv_pid_pointer (p, pid));
}

static void rcv_pid_append_pids (RcvPid * p, Array * pids)
{

  for (int i = 0; i < p->npid; i++) {
    int pid = p->rcv[i].pid, j, * a;
    for (j = 0, a = pids->p; j < pids->len/sizeof(int); j++,a++)
      if (*a == pid)
 break;
    if (j == pids->len/sizeof(int))
      array_append (pids, &pid, sizeof(int));
  }
}

void rcv_pid_write (RcvPid * p, const char * name)
{
  for (int i = 0; i < p->npid; i++) {
    Rcv * rcv = &p->rcv[i];
    char fname[80];
    sprintf (fname, "%s-%d-%d", name, pid(), rcv->pid);
    FILE * fp = fopen (fname, "w");
    rcv_print (rcv, fp, "");
    fclose (fp);
  }
}

static void rcv_pid_print (RcvPid * p, FILE * fp, const char * prefix)
{
  for (int i = 0; i < p->npid; i++)
    rcv_print (&p->rcv[i], fp, prefix);
}

static void rcv_pid_destroy (RcvPid * p)
{
  for (int i = 0; i < p->npid; i++)
    rcv_destroy (&p->rcv[i]);
  pfree (p->rcv,__func__,__FILE__,__LINE__);
  pfree (p->name,__func__,__FILE__,__LINE__);
  pfree (p,__func__,__FILE__,__LINE__);
}

static Boundary * mpi_boundary = NULL;






void debug_mpi (FILE * fp1);

static void apply_bc (Rcv * rcv, scalar * list, vector * listf, int l,
        MPI_Status s)
{
  double * b = rcv->buf;
   { foreach_cache_level(rcv->halo[l], l,){

#line 165 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"
 {
    if (list) for (scalar s = *list, *_i58 = list; ((scalar *)&s)->i >= 0; s = *++_i58)
      val(s,0,0,0) = *b++;
    if (listf) for (vector v = *listf, *_i59 = listf; ((scalar *)&v)->i >= 0; v = *++_i59) {
      {
#line 169
 {
 val(v.x,0,0,0) = *b++;
 if (allocated(1,0,0))
   val(v.x,1,0,0) = *b;
 b++;
      }
#line 169
 {
 val(v.y,0,0,0) = *b++;
 if (allocated(0,1,0))
   val(v.y,0,1,0) = *b;
 b++;
      }
#line 169
 {
 val(v.z,0,0,0) = *b++;
 if (allocated(0,0,1))
   val(v.z,0,0,1) = *b;
 b++;
      }}
    }
  } } end_foreach_cache_level(); }
  size_t size = b - (double *) rcv->buf;
  pfree (rcv->buf,__func__,__FILE__,__LINE__);
  rcv->buf = NULL;

  int rlen;
  MPI_Get_count (&s, MPI_DOUBLE, &rlen);
  if (rlen != size) {
    fprintf (qstderr(),
      "rlen (%d) != size (%ld), %d receiving from %d at level %d\n"
      "Calling debug_mpi(NULL)...\n"
      "Aborting...\n",
      rlen, size, pid(), rcv->pid, l);
    fflush (qstderr());
    debug_mpi (NULL);
    MPI_Abort (MPI_COMM_WORLD, -2);
  }
}
#line 215 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"
static void mpi_recv_check (void * buf, int count, MPI_Datatype datatype,
       int source, int tag,
       MPI_Comm comm, MPI_Status * status,
       const char * name)
{
#line 250 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"
  int errorcode = MPI_Recv (buf, count, datatype, source, tag, comm, status);
  if (errorcode != MPI_SUCCESS) {
    char string[MPI_MAX_ERROR_STRING];
    int resultlen;
    MPI_Error_string (errorcode, string, &resultlen);
    fprintf (qstderr(),
      "ERROR MPI_Recv \"%s\" (count = %d, source = %d, tag = %d):\n%s\n"
      "Calling debug_mpi(NULL)...\n"
      "Aborting...\n",
      name, count, source, tag, string);
    fflush (qstderr());
    debug_mpi (NULL);
    MPI_Abort (MPI_COMM_WORLD, -1);
  }





}

static void rcv_pid_receive (RcvPid * m, scalar * list, vector * listf, int l)
{
  if (m->npid == 0)
    return;

  prof_start ("rcv_pid_receive");

  int len = list_len (list) + 2*3*vectors_len (listf);

  MPI_Request r[m->npid];
  Rcv * rrcv[m->npid];
  int nr = 0;
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      assert (!rcv->buf);
      rcv->buf = pmalloc (sizeof (double)*rcv->halo[l].n*len,__func__,__FILE__,__LINE__);






      MPI_Irecv (rcv->buf, rcv->halo[l].n*len, MPI_DOUBLE, rcv->pid,
   (l), MPI_COMM_WORLD, &r[nr]);
      rrcv[nr++] = rcv;






    }
  }


  if (nr > 0) {
    int i;
    MPI_Status s;
    MPI_Waitany (nr, r, &i, &s);
    while (i != MPI_UNDEFINED) {
      Rcv * rcv = rrcv[i];
      assert (l <= rcv->depth && rcv->halo[l].n > 0);
      assert (rcv->buf);
      apply_bc (rcv, list, listf, l, s);
      MPI_Waitany (nr, r, &i, &s);
    }
  }

  prof_stop();
}

static void rcv_pid_wait (RcvPid * m)
{

  for (int i = 0; i < m->npid; i++)
    rcv_free_buf (&m->rcv[i]);
}

static void rcv_pid_send (RcvPid * m, scalar * list, vector * listf, int l)
{
  if (m->npid == 0)
    return;

  prof_start ("rcv_pid_send");

  int len = list_len (list) + 2*3*vectors_len (listf);


  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      assert (!rcv->buf);
      rcv->buf = pmalloc (sizeof (double)*rcv->halo[l].n*len,__func__,__FILE__,__LINE__);
      double * b = rcv->buf;
       { foreach_cache_level(rcv->halo[l], l,){

#line 346 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"
 {
 if (list) for (scalar s = *list, *_i60 = list; ((scalar *)&s)->i >= 0; s = *++_i60)
   *b++ = val(s,0,0,0);
 if (listf) for (vector v = *listf, *_i61 = listf; ((scalar *)&v)->i >= 0; v = *++_i61)
   {
#line 350
 {
     *b++ = val(v.x,0,0,0);
     *b++ = allocated(1,0,0) ? val(v.x,1,0,0) : undefined;
   }
#line 350
 {
     *b++ = val(v.y,0,0,0);
     *b++ = allocated(0,1,0) ? val(v.y,0,1,0) : undefined;
   }
#line 350
 {
     *b++ = val(v.z,0,0,0);
     *b++ = allocated(0,0,1) ? val(v.z,0,0,1) : undefined;
   }}
      } } end_foreach_cache_level(); }





      MPI_Isend (rcv->buf, (b - (double *) rcv->buf),
   MPI_DOUBLE, rcv->pid, (l), MPI_COMM_WORLD,
   &rcv->r);
    }
  }

  prof_stop();
}

static void rcv_pid_sync (SndRcv * m, scalar * list, int l)
{
  scalar * listr = NULL;
  vector * listf = NULL;
  if (list) for (scalar s = *list, *_i62 = list; ((scalar *)&s)->i >= 0; s = *++_i62)
    if (!is_constant(s)) {
      if (_attribute[s.i].face)
 listf = vectors_add (listf, _attribute[s.i].v);
      else
 listr = list_add (listr, s);
    }
  rcv_pid_send (m->snd, listr, listf, l);
  rcv_pid_receive (m->rcv, listr, listf, l);
  rcv_pid_wait (m->snd);
  pfree (listr,__func__,__FILE__,__LINE__);
  pfree (listf,__func__,__FILE__,__LINE__);
}

static void snd_rcv_destroy (SndRcv * m)
{
  rcv_pid_destroy (m->rcv);
  rcv_pid_destroy (m->snd);
}

static void snd_rcv_init (SndRcv * m, const char * name)
{
  char s[strlen(name) + 5];
  strcpy (s, name);
  strcat (s, ".rcv");
  m->rcv = rcv_pid_new (s);
  strcpy (s, name);
  strcat (s, ".snd");
  m->snd = rcv_pid_new (s);
}

static void mpi_boundary_destroy (Boundary * b)
{
  MpiBoundary * m = (MpiBoundary *) b;
  snd_rcv_destroy (&m->restriction);
  snd_rcv_destroy (&m->restriction_root);
  snd_rcv_destroy (&m->halo_restriction);
  snd_rcv_destroy (&m->prolongation);
  array_free (m->send);
  array_free (m->receive);
  pfree (m,__func__,__FILE__,__LINE__);
}


static void mpi_boundary_restriction (const Boundary * b, scalar * list, int l)
{ trace ("mpi_boundary_restriction", "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h", 418);
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->restriction, list, l);
  rcv_pid_sync (&m->restriction_root, list, l);
 end_trace("mpi_boundary_restriction", "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h", 422); }


static void mpi_boundary_halo_restriction (const Boundary * b,
        scalar * list, int l)
{ trace ("mpi_boundary_halo_restriction", "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h", 427);
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->halo_restriction, list, l);
 end_trace("mpi_boundary_halo_restriction", "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h", 430); }


static void mpi_boundary_halo_prolongation (const Boundary * b,
         scalar * list, int l, int depth)
{ trace ("mpi_boundary_halo_prolongation", "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h", 435);
  MpiBoundary * m = (MpiBoundary *) b;
  if (true ) {
    rcv_pid_sync (&m->restriction, list, l);
    rcv_pid_sync (&m->restriction_root, list, l);
  }
  else
    rcv_pid_sync (&m->prolongation, list, l);
 end_trace("mpi_boundary_halo_prolongation", "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h", 443); }

void mpi_boundary_new()
{
  mpi_boundary = pcalloc (1, sizeof (MpiBoundary),__func__,__FILE__,__LINE__);
  mpi_boundary->destroy = mpi_boundary_destroy;
  mpi_boundary->restriction = mpi_boundary_restriction;
  mpi_boundary->halo_restriction = mpi_boundary_halo_restriction;
  mpi_boundary->halo_prolongation = mpi_boundary_halo_prolongation;
  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;
  snd_rcv_init (&mpi->restriction, "restriction");
  snd_rcv_init (&mpi->restriction_root, "restriction_root");
  snd_rcv_init (&mpi->halo_restriction, "halo_restriction");
  snd_rcv_init (&mpi->prolongation, "prolongation");
  mpi->send = array_new();
  mpi->receive = array_new();
  add_boundary (mpi_boundary);
}

static FILE * fopen_prefix (FILE * fp, const char * name, char * prefix)
{
  if (fp) {
    sprintf (prefix, "%s-%d ", name, pid());
    return fp;
  }
  else {
    strcpy (prefix, "");
    char fname[80];
    if (debug_iteration >= 0)
      sprintf (fname, "%s-%d-%d", name, debug_iteration, pid());
    else
      sprintf (fname, "%s-%d", name, pid());
    return fopen (fname, "w");
  }
}

void debug_mpi (FILE * fp1)
{
  void output_cells (FILE * fp);

  char prefix[80];
  FILE * fp;


  if (fp1 == NULL) {
    char name[80];
    sprintf (name, "halo-%d", pid()); remove (name);
    sprintf (name, "cells-%d", pid()); remove (name);
    sprintf (name, "faces-%d", pid()); remove (name);
    sprintf (name, "neighbors-%d", pid()); remove (name);
    sprintf (name, "restriction-%d", pid()); remove (name);
    sprintf (name, "restriction-root-%d", pid()); remove (name);
    sprintf (name, "mpi-restriction-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-halo-restriction-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-prolongation-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-restriction-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-halo-restriction-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-prolongation-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-border-%d", pid()); remove (name);
    sprintf (name, "exterior-%d", pid()); remove (name);
    sprintf (name, "depth-%d", pid()); remove (name);
    sprintf (name, "refined-%d", pid()); remove (name);
  }


  fp = fopen_prefix (fp1, "halo", prefix);
  for (int l = 0; l < depth(); l++)
     { foreach_halo (prolongation, l){

#line 510 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"

       { foreach_child()
      fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); end_foreach_child(); }; } end_foreach_halo(); }
  if (!fp1)
    fclose (fp);

  if (!fp1) {
    fp = fopen_prefix (fp1, "cells", prefix);
    output_cells (fp);
    fclose (fp);
  }

  fp = fopen_prefix (fp1, "faces", prefix);
   { foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 523
{

#line 523 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"

    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 523
{

#line 523 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"

    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 523
{

#line 523 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"

    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); }  }}  end_foreach_face_generic()
#line 524
 end_foreach_face(); }
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "neighbors", prefix);
   { foreach(){

#line 529 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"
 {
    int n = 0;
     { foreach_neighbor(1)
      if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
 n++; end_foreach_neighbor(); }
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, cell.neighbors);
    assert (cell.neighbors == n);
  } } end_foreach(); }
  if (!fp1)
    fclose (fp);


  fp = fopen_prefix (fp1, "restriction", prefix);
  for (int l = 0; l < depth(); l++)
     { foreach_halo (restriction, l){

#line 543 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"

      fprintf (fp, "%s%g %g %g %d %d\n", prefix, x, y, z, level, cell.neighbors); } end_foreach_halo(); }
  if (!fp1)
    fclose (fp);

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;

  fp = fopen_prefix (fp1, "mpi-restriction-rcv", prefix);
  rcv_pid_print (mpi->restriction.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-restriction-root-rcv", prefix);
  rcv_pid_print (mpi->restriction_root.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-halo-restriction-rcv", prefix);
  rcv_pid_print (mpi->halo_restriction.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-prolongation-rcv", prefix);
  rcv_pid_print (mpi->prolongation.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-restriction-snd", prefix);
  rcv_pid_print (mpi->restriction.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-restriction-root-snd", prefix);
  rcv_pid_print (mpi->restriction_root.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-halo-restriction-snd", prefix);
  rcv_pid_print (mpi->halo_restriction.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-prolongation-snd", prefix);
  rcv_pid_print (mpi->prolongation.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-border", prefix);
   { foreach_cell(){

#line 591 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"
 {
    if (is_border(cell))
      fprintf (fp, "%s%g %g %g %d %d %d\n",
        prefix, x, y, z, level, cell.neighbors, cell.pid);
    else
      continue;
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "exterior", prefix);
   { foreach_cell(){

#line 604 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"
 {
    if (!is_local(cell))
      fprintf (fp, "%s%g %g %g %d %d %d %d\n",
        prefix, x, y, z, level, cell.neighbors,
        cell.pid, cell.flags & leaf);






  } } end_foreach_cell(); }
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "depth", prefix);
  fprintf (fp, "depth: %d %d\n", pid(), depth());
  fprintf (fp, "======= restriction.snd ======\n");
  RcvPid * snd = mpi->restriction.snd;
  for (int i = 0; i < snd->npid; i++)
    fprintf (fp, "%d %d %d\n", pid(), snd->rcv[i].pid, snd->rcv[i].maxdepth);
  fprintf (fp, "======= restriction.rcv ======\n");
  snd = mpi->restriction.rcv;
  for (int i = 0; i < snd->npid; i++)
    fprintf (fp, "%d %d %d\n", pid(), snd->rcv[i].pid, snd->rcv[i].maxdepth);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "refined", prefix);
   { foreach_cache (((Quadtree *)grid)->refined,){

#line 633 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"

    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); } end_foreach_cache(); }
  if (!fp1)
    fclose (fp);
}

static void snd_rcv_free (SndRcv * p)
{
  char name[strlen(p->rcv->name) + 1];
  strcpy (name, p->rcv->name);
  rcv_pid_destroy (p->rcv);
  p->rcv = rcv_pid_new (name);
  strcpy (name, p->snd->name);
  rcv_pid_destroy (p->snd);
  p->snd = rcv_pid_new (name);
}

static bool is_root (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 651 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"

  if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
     { foreach_child()
      if (is_local(cell))
 return true; end_foreach_child(); }
  return false;
}

static void append_pid (Array * pids, Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 660 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"

  for (int i = 0, * p = (int *) pids->p; i < pids->len/sizeof(int); i++, p++)
    if (*p == cell.pid)
      return;
  array_append (pids, &cell.pid, sizeof(int));
}


static bool is_local_prolongation (Point point, Point p)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 669 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"




  struct { int x, y, z; } dp = {p.i - point.i, p.j - point.j, p.k - point.k};

  {
#line 675
 {
    if (dp.x == 0 && ((!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) || (!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0)))
      return true;
    if ((!is_leaf (neighbor(dp.x,0,0)) && neighbor(dp.x,0,0).neighbors && neighbor(dp.x,0,0).pid >= 0))
      return true;
  }
#line 675
 {
    if (dp.y == 0 && ((!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) || (!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0)))
      return true;
    if ((!is_leaf (neighbor(0,dp.y,0)) && neighbor(0,dp.y,0).neighbors && neighbor(0,dp.y,0).pid >= 0))
      return true;
  }
#line 675
 {
    if (dp.z == 0 && ((!is_leaf (neighbor(0,0,-1)) && neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0) || (!is_leaf (neighbor(0,0,1)) && neighbor(0,0,1).neighbors && neighbor(0,0,1).pid >= 0)))
      return true;
    if ((!is_leaf (neighbor(0,0,dp.z)) && neighbor(0,0,dp.z).neighbors && neighbor(0,0,dp.z).pid >= 0))
      return true;
  }}
  return false;
}



static int locals_pids (Point point, Array * pids)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 687 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"

  if (is_leaf(cell)) {
    if (is_local(cell)) {
      Point p = point;
       { foreach_neighbor(1) {
 if ((cell.pid >= 0 && cell.pid != pid()) &&
     ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) || is_local_prolongation (point, p)))
   append_pid (pids, point);
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
    { foreach_child()
     if ((cell.pid >= 0 && cell.pid != pid()))
       append_pid (pids, point); end_foreach_child(); }
      } end_foreach_neighbor(); }
    }
  }
  else
     { foreach_neighbor(1) {
      if ((cell.pid >= 0 && cell.pid != pid()))
 append_pid (pids, point);
      if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
  { foreach_child()
   if ((cell.pid >= 0 && cell.pid != pid()))
     append_pid (pids, point); end_foreach_child(); }
    } end_foreach_neighbor(); }
  return pids->len/sizeof(int);
}

static int root_pids (Point point, Array * pids)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 715 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"

   { foreach_child()
    if ((cell.pid >= 0 && cell.pid != pid()))
      append_pid (pids, point); end_foreach_child(); }
  return pids->len/sizeof(int);
}







static void rcv_pid_row (RcvPid * m, int l, int * row)
{
  for (int i = 0; i < npe(); i++)
    row[i] = 0;
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0)
      row[rcv->pid] = rcv->halo[l].n;
  }
}

void check_snd_rcv_matrix (SndRcv * sndrcv, const char * name)
{
  int maxlevel = depth();
  mpi_all_reduce (maxlevel, MPI_INT, MPI_MAX);
  int * row = pmalloc (npe()*sizeof(int),__func__,__FILE__,__LINE__);
  for (int l = 0; l <= maxlevel; l++) {
    int status = 0;
    if (pid() == 0) {


      int ** send = matrix_new (npe(), npe(), sizeof(int));
      int ** receive = matrix_new (npe(), npe(), sizeof(int));
      rcv_pid_row (sndrcv->snd, l, row);
      MPI_Gather (row, npe(), MPI_INT, &send[0][0], npe(), MPI_INT, 0,
    MPI_COMM_WORLD);
      rcv_pid_row (sndrcv->rcv, l, row);
      MPI_Gather (row, npe(), MPI_INT, &receive[0][0], npe(), MPI_INT, 0,
    MPI_COMM_WORLD);

      int * astatus = pmalloc (npe()*sizeof(int),__func__,__FILE__,__LINE__);
      for (int i = 0; i < npe(); i++)
 astatus[i] = 0;
      for (int i = 0; i < npe(); i++)
 for (int j = 0; j < npe(); j++)
   if (send[i][j] != receive[j][i]) {
     fprintf (qstderr(), "%s: %d sends    %d to   %d at level %d\n",
       name, i, send[i][j], j, l);
     fprintf (qstderr(), "%s: %d receives %d from %d at level %d\n",
       name, j, receive[j][i], i, l);
     fflush (qstderr());
     for (int k = i - 2; k <= i + 2; k++)
       if (k >= 0 && k < npe())
  astatus[k] = 1;
     for (int k = j - 2; k <= j + 2; k++)
       if (k >= 0 && k < npe())
  astatus[k] = 1;
   }
      MPI_Scatter (astatus, 1, MPI_INT, &status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      pfree (astatus,__func__,__FILE__,__LINE__);

      matrix_free (send);
      matrix_free (receive);
    }
    else {
      rcv_pid_row (sndrcv->snd, l, row);
      MPI_Gather (row, npe(), MPI_INT, NULL, npe(), MPI_INT, 0, MPI_COMM_WORLD);
      rcv_pid_row (sndrcv->rcv, l, row);
      MPI_Gather (row, npe(), MPI_INT, NULL, npe(), MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Scatter (NULL, 1, MPI_INT, &status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    if (status) {
      fprintf (qstderr(),
        "check_snd_rcv_matrix \"%s\" failed\n"
        "Calling debug_mpi(NULL)...\n"
        "Aborting...\n",
        name);
      fflush (qstderr());
      debug_mpi (NULL);
      MPI_Abort (MPI_COMM_WORLD, -3);
    }
  }
  pfree (row,__func__,__FILE__,__LINE__);
}


void mpi_boundary_update()
{ trace ("mpi_boundary_update", "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h", 805);
  if (npe() == 1)
    { ; end_trace("mpi_boundary_update", "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h", 807);  return; }

  prof_start ("mpi_boundary_update");

  MpiBoundary * m = (MpiBoundary *) mpi_boundary;
  SndRcv * restriction = &m->restriction;
  SndRcv * restriction_root = &m->restriction_root;
  SndRcv * halo_restriction = &m->halo_restriction;

  snd_rcv_free (restriction);
  snd_rcv_free (restriction_root);
  snd_rcv_free (halo_restriction);

  static const unsigned short used = 1 << user;
   { foreach_cell(){

#line 821 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"
 {
    if (is_active(cell) && !is_border(cell))



      continue;
    if (cell.neighbors) {

      Array pids = {NULL, 0, 0};
      int n = locals_pids (point, &pids);
      if (n) {
  { foreach_child()
   if (is_local(cell))
     for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
       rcv_pid_append (restriction->snd, *p, point); end_foreach_child(); }
 pfree (pids.p,__func__,__FILE__,__LINE__);
      }

      bool locals = false;
      if (is_leaf(cell)) {
 if ((cell.pid >= 0 && cell.pid != pid())) {
   Point p = point;
    { foreach_neighbor(1)
     if ((is_local(cell) &&
   ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) || is_local_prolongation (point, p))) ||
  is_root(point))
       locals = true, foreach_neighbor_break(); end_foreach_neighbor(); }
 }
      }
      else
  { foreach_neighbor(1)
   if (is_local(cell) || is_root(point))
     locals = true, foreach_neighbor_break(); end_foreach_neighbor(); }
      if (locals)
  { foreach_child()
   if ((cell.pid >= 0 && cell.pid != pid()))
            rcv_pid_append (restriction->rcv, cell.pid, point),
       cell.flags |= used; end_foreach_child(); }


      if (!is_leaf(cell)) {

 if (is_local(cell)) {
   Array pids = {NULL, 0, 0};

   int n = root_pids (point, &pids);
   if (n) {
      { foreach_neighbor()
       for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
  if (cell.pid >= 0 && cell.pid != *p)
    rcv_pid_append (restriction_root->snd, *p, point); end_foreach_neighbor(); }
     pfree (pids.p,__func__,__FILE__,__LINE__);
   }
 }

 else if ((cell.pid >= 0 && cell.pid != pid())) {
   bool root = false;
    { foreach_child()
     if (is_local(cell))
       root = true, foreach_child_break(); end_foreach_child(); }
   if (root) {
     int pid = cell.pid;
      { foreach_neighbor()
       if ((cell.pid >= 0 && cell.pid != pid()))
  rcv_pid_append (restriction_root->rcv, pid, point),
    cell.flags |= used; end_foreach_neighbor(); }
   }
 }
      }
    }
  } } end_foreach_cell(); }





  static const unsigned short keep = 1 << (user + 1);
  for (int l = depth(); l >= 0; l--)
     { foreach_cell(){

#line 899 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"

      if (level == l) {
 if (level > 0 && (cell.pid < 0 || is_local(cell) || (cell.flags & used)))
   aparent(0,0,0).flags |= keep;
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) && !(cell.flags & keep))
   coarsen_cell (point, NULL);
 cell.flags &= ~(used|keep);
 continue;
      } } end_foreach_cell(); }


  m->send->len = m->receive->len = 0;
  rcv_pid_append_pids (restriction->snd, m->send);
  rcv_pid_append_pids (restriction_root->snd, m->send);
  rcv_pid_append_pids (restriction->rcv, m->receive);
  rcv_pid_append_pids (restriction_root->rcv, m->receive);

  prof_stop();
#line 932 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"
  scalar halov= new_scalar("halov");
  ((Quadtree *)grid)->dirty = true;
  for (int l = 0; l <= depth(); l++) {
     { foreach_cell(){

#line 935 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"
 {
      if (level == l) {
 if (is_local(cell)) {
   bool restriction = level > 0 && coarse(halov,0,0,0);
   if (restriction && !is_local(aparent(0,0,0)))
     rcv_pid_append (halo_restriction->snd, aparent(0,0,0).pid, point);
   if (!is_leaf(cell)) {
     if (!restriction)
        { foreach_neighbor()
  if (is_leaf(cell) && !(cell.pid < 0))
    restriction = true, foreach_neighbor_break(); end_foreach_neighbor(); }
     if (restriction) {
       cell.flags |= halo;
       val(halov,0,0,0) = true;
     }
     else {
       cell.flags &= ~halo;
       val(halov,0,0,0) = false;
     }
   }
 }
 else {
   bool restriction = level > 0 && (aparent(0,0,0).flags & halo);
   if (restriction && is_local(aparent(0,0,0)))
     rcv_pid_append (halo_restriction->rcv, cell.pid, point);
   if (!is_leaf(cell)) {
     if (!restriction)

        { foreach_neighbor()
  if (allocated(0,0,0) && is_leaf(cell) && is_local(cell))
    restriction = true, foreach_neighbor_break(); end_foreach_neighbor(); }
     if (restriction)
       cell.flags |= halo;
     else
       cell.flags &= ~halo;
   }
 }
 continue;
      }
      if (is_leaf(cell))
 continue;
    } } end_foreach_cell(); }

    mpi_boundary_restriction (mpi_boundary, ((scalar []){halov,{-1}}), l);
  }
 delete (((scalar []){halov,{-1}}));  end_trace("mpi_boundary_update", "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h", 980); }


void mpi_boundary_refine (scalar * list, int file)
{ trace ("mpi_boundary_refine", "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h", 984);
  prof_start ("mpi_boundary_refine");

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;


  Array * snd = mpi->send;
  MPI_Request r[2*snd->len/sizeof(int)];
  int nr = 0;
  for (int i = 0, * dest = snd->p; i < snd->len/sizeof(int); i++,dest++) {
    int len = ((Quadtree *)grid)->refined.n;
    MPI_Isend (&((Quadtree *)grid)->refined.n, 1, MPI_INT, *dest,
        (128), MPI_COMM_WORLD, &r[nr++]);
    if (len > 0)
      MPI_Isend (((Quadtree *)grid)->refined.p, sizeof(Index)/sizeof(int)*len,
   MPI_INT, *dest, (128), MPI_COMM_WORLD, &r[nr++]);
  }



  Array * rcv = mpi->receive;
  Cache rerefined = {NULL, 0, 0};
  for (int i = 0, * source = rcv->p; i < rcv->len/sizeof(int); i++,source++) {
    int len;
    mpi_recv_check (&len, 1, MPI_INT, *source, (128),
      MPI_COMM_WORLD, MPI_STATUS_IGNORE,
      "mpi_boundary_refine (len)");
    if (len > 0) {
      Index p[len];
      mpi_recv_check (p, sizeof(Index)/sizeof(int)*len,
        MPI_INT, *source, (128),
        MPI_COMM_WORLD, MPI_STATUS_IGNORE,
        "mpi_boundary_refine (p)");
      Cache refined = {p, len, len};
       { foreach_cache (refined,){

#line 1018 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"

 if (level <= depth() && allocated(0,0,0)) {
   if (is_leaf(cell)) {
     bool neighbors = false;
      { foreach_neighbor()
       if (allocated(0,0,0) && (is_active(cell) || is_local(aparent(0,0,0))))
  neighbors = true, foreach_neighbor_break(); end_foreach_neighbor(); }

     if (neighbors)
       refine_cell (point, list, 0, &rerefined);
   }
 } } end_foreach_cache(); }
    }
  }


  if (nr)
    MPI_Waitall (nr, r, MPI_STATUSES_IGNORE);


  pfree (((Quadtree *)grid)->refined.p,__func__,__FILE__,__LINE__);
  ((Quadtree *)grid)->refined = rerefined;

  prof_stop();



  mpi_all_reduce (rerefined.n, MPI_INT, MPI_SUM);
  if (rerefined.n)
    mpi_boundary_refine (list, 0);
 end_trace("mpi_boundary_refine", "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h", 1048); }

static void check_depth()
{
#line 1083 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"
}

typedef struct {
  int refined, leaf;
} Remote;




void mpi_boundary_coarsen (int l, int too_fine)
{ trace ("mpi_boundary_coarsen", "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h", 1093);
  if (npe() == 1)
    { ; end_trace("mpi_boundary_coarsen", "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h", 1095);  return; }

  check_depth();

  assert (sizeof(Remote) == sizeof(double));

  scalar remote= new_scalar("remote");
   { foreach_cell(){

#line 1102 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"
 {
    if (level == l) {
      if (is_local(cell)) {
 ((Remote *)&val(remote,0,0,0))->refined = (!is_leaf (cell) && cell.neighbors && cell.pid >= 0);
 ((Remote *)&val(remote,0,0,0))->leaf = is_leaf(cell);
      }
      else {
 ((Remote *)&val(remote,0,0,0))->refined = true;
 ((Remote *)&val(remote,0,0,0))->leaf = false;
      }
      continue;
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }
  mpi_boundary_restriction (mpi_boundary, ((scalar []){remote,{-1}}), l);

   { foreach_cell(){

#line 1119 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"
 {
    if (level == l) {
      if (!is_local(cell)) {
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) && !((Remote *)&val(remote,0,0,0))->refined)
   coarsen_cell_recursive (point, NULL);
 else if (is_leaf(cell) && cell.neighbors && ((Remote *)&val(remote,0,0,0))->leaf) {
   int pid = cell.pid;
    { foreach_child()
     cell.pid = pid; end_foreach_child(); }
 }
      }
      continue;
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  check_depth();

  if (l > 0) {
     { foreach_cell(){

#line 1139 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"
 {
      if (level == l) {
 val(remote,0,0,0) = is_local(cell) ? cell.neighbors : 0;
 continue;
      }
      if (is_leaf(cell))
 continue;
    } } end_foreach_cell(); }
    mpi_boundary_restriction (mpi_boundary, ((scalar []){remote,{-1}}), l);
     { foreach_cell(){

#line 1148 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"
 {
      if (level == l)
 if (!is_local(cell) && is_local(aparent(0,0,0)) && val(remote,0,0,0)) {
   aparent(0,0,0).flags &= ~too_fine;
   continue;
 }
      if (is_leaf(cell))
 continue;
    } } end_foreach_cell(); }
  }
 delete (((scalar []){remote,{-1}}));  end_trace("mpi_boundary_coarsen", "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h", 1158); }

static void flag_border_cells()
{
   { foreach_cell(){

#line 1162 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"
 {
    if (is_active(cell)) {
      short flags = cell.flags & ~border;
       { foreach_neighbor() {
 if (!is_local(cell) || (level > 0 && !is_local(aparent(0,0,0))))
   flags |= border, foreach_neighbor_break();

 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
    { foreach_child()
     if (!is_local(cell))
       flags |= border, foreach_child_break(); end_foreach_child(); }
 if (flags & border)
   foreach_neighbor_break();
      } end_foreach_neighbor(); }
      cell.flags = flags;
    }
    else {
      cell.flags &= ~border;

    }
    if (is_leaf(cell)) {
      if (cell.neighbors) {
  { foreach_child()
   cell.flags &= ~border; end_foreach_child(); }
 if (is_border(cell)) {
   bool remote = false;
    { foreach_neighbor (2/2)
     if (!is_local(cell))
       remote = true, foreach_neighbor_break(); end_foreach_neighbor(); }
   if (remote)
      { foreach_child()
       cell.flags |= border; end_foreach_child(); }
 }
      }
      continue;
    }
  } } end_foreach_cell(); }
}

static int balanced_pid (long index, long nt, int nproc)
{
  long ne = max(1, nt/nproc), nr = nt % nproc;
  int pid = index < nr*(ne + 1) ?
    index/(ne + 1) :
    nr + (index - nr*(ne + 1))/ne;
  return min(nproc - 1, pid);
}


void mpi_partitioning()
{ trace ("mpi_partitioning", "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h", 1212);
  prof_start ("mpi_partitioning");

  long nt = 0;
   { foreach(){

#line 1216 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"

    nt++; } end_foreach(); }


  long i = 0;
  ((Quadtree *)grid)->dirty = true;
   { foreach_cell_post (is_active (cell)){

#line 1222 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"

    if (is_active (cell)) {
      if (is_leaf (cell)) {
 cell.pid = balanced_pid (i++, nt, npe());
 if (cell.neighbors > 0) {
   int pid = cell.pid;
    { foreach_child()
     cell.pid = pid; end_foreach_child(); }
 }
 if (!is_local(cell))
   cell.flags &= ~active;
      }
      else {
 cell.pid = child(0,0,0).pid;
 bool inactive = true;
  { foreach_child()
   if (is_active(cell))
     inactive = false, foreach_child_break(); end_foreach_child(); }
 if (inactive)
   cell.flags &= ~active;
      }
    } } end_foreach_cell_post(); }

  flag_border_cells();

  prof_stop();

  mpi_boundary_update();
 end_trace("mpi_partitioning", "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h", 1250); }
#line 1272 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"

double z_indexing (scalar index, bool leaves)
{ trace ("z_indexing", "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h", 1274);






  scalar size= new_scalar("size");




   { foreach(){

#line 1286 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"

    val(size,0,0,0) = 1; } end_foreach(); }






  { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){size,{-1}}), depth()); };
  for (int l = depth() - 1; l >= 0; l--) {
     { foreach_coarse_level(l){

#line 1296 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"
 {
      double sum = !leaves;
       { foreach_child()
 sum += val(size,0,0,0); end_foreach_child(); }
      val(size,0,0,0) = sum;
    } } end_foreach_coarse_level(); }
    { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){size,{-1}}), l); };
  }






  double maxi = -1.;
  if (pid() == 0)
     { foreach_level(0){

#line 1312 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"

      maxi = val(size,0,0,0) - 1.; } end_foreach_level(); }




   { foreach_level(0){

#line 1318 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"

    val(index,0,0,0) = 0; } end_foreach_level(); }
  for (int l = 0; l < depth(); l++) {

    { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){index,{-1}}), l); };
     { foreach_cell(){

#line 1323 "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h"
 {
      if (level == l) {
 if (is_leaf(cell)) {
   if (is_local(cell) && cell.neighbors) {
     int i = val(index,0,0,0);
      { foreach_child()
       val(index,0,0,0) = i; end_foreach_child(); }
   }
 }
 else {
   bool loc = is_local(cell);
   if (!loc)
      { foreach_child()
       if (is_local(cell))
  loc = true, foreach_child_break(); end_foreach_child(); }
   if (loc) {
     int i = val(index,0,0,0) + !leaves;
      { foreach_child() {
       val(index,0,0,0) = i;
       i += val(size,0,0,0);
     } end_foreach_child(); }
   }
 }
 continue;
      }
      if (is_leaf(cell))
 continue;
    } } end_foreach_cell(); }
  }

  { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){index,{-1}}), depth()); };

  { double _ret =  maxi; delete (((scalar []){size,{-1}}));  end_trace("z_indexing", "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h", 1355);  return _ret; }
 delete (((scalar []){size,{-1}}));  end_trace("z_indexing", "/home/popinet/basilisk-octree/src/grid/quadtree-mpi.h", 1356); }
#line 1716 "/home/popinet/basilisk-octree/src/grid/tree.h"
#line 1 "grid/balance.h"
#line 1 "/home/popinet/basilisk-octree/src/grid/balance.h"


typedef struct {
  short leaf, prolongation;
  int pid;
} NewPid;



#if TRASH
# define is_newpid() (!isnan(val(newpid,0,0,0)) && ((NewPid *)&val(newpid,0,0,0))->pid > 0)
#else
# define is_newpid() (((NewPid *)&val(newpid,0,0,0))->pid > 0)
#endif

Array * tree (size_t size, scalar newpid)
{
  const unsigned short sent = 1 << user, next = 1 << (user + 1);
  Array * a = array_new();

   { foreach_cell_post_all (true){

#line 21 "/home/popinet/basilisk-octree/src/grid/balance.h"

    if (level > 0 && (cell.flags & (sent|next)))
      aparent(0,0,0).flags |= next; } end_foreach_cell_post_all(); }

  bool empty = true;
   { foreach_cell_all(){

#line 26 "/home/popinet/basilisk-octree/src/grid/balance.h"
 {
    if (cell.flags & sent) {
      array_append (a, &cell, size);
      cell.flags &= ~sent;
      empty = false;
    }
    else {
      if (cell.pid >= 0 && ((NewPid *)&val(newpid,0,0,0))->leaf)
 assert (is_leaf(cell));
      if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0)) {


 bool prolo = false;
  { foreach_child()
   if (((NewPid *)&val(newpid,0,0,0))->prolongation)
     prolo = true; end_foreach_child(); }
 if (prolo) {

   cell.flags |= leaf;
   array_append (a, &cell, sizeof(Cell));
   cell.flags &= ~leaf;
 }
 else
   array_append (a, &cell, sizeof(Cell));
      }
      else
 array_append (a, &cell, sizeof(Cell));
    }
    if (cell.flags & next)
      cell.flags &= ~next;
    else
      continue;
  } } end_foreach_cell_all(); }

  if (empty)
    a->len = 0;
  return a;
}

#define foreach_tree(t, size, list)\
{\
  const unsigned short _sent = 1 << user, _next = 1 << (user + 1);\
  scalar * _list = list;\
  char * _i = (char *) (t)->p;\
  foreach_cell_all() {\
    Cell * c = (Cell *) _i;\
    if (c->flags & _sent) {\
      _i += size;\

#line 74


#define end_foreach_tree()\
    }\
    else\
      _i += sizeof(Cell);\
    if (c->flags & _next) {\
      assert (c->neighbors);\
      if (!(c->flags & leaf) && is_leaf(cell) &&\
   (!is_newpid() || !((NewPid *)&val(newpid,0,0,0))->leaf))\
\
 refine_cell (point, _list, 0, NULL);\
      else if (!cell.neighbors)\
\
 alloc_children (point);\
    }\
    else\
      continue;\
  } end_foreach_cell_all();\
}\

#line 94


Array * neighborhood (scalar newpid, int nextpid, FILE * fp)
{
  const unsigned short sent = 1 << user;
   { foreach_cell(){

#line 99 "/home/popinet/basilisk-octree/src/grid/balance.h"
 {

    bool root = false;
    if ((!is_local(cell) || ((NewPid *)&val(newpid,0,0,0))->pid - 1 != nextpid) && (!is_leaf (cell) && cell.neighbors && cell.pid >= 0)) {
       { foreach_child()
 if (is_local(cell) && ((NewPid *)&val(newpid,0,0,0))->pid - 1 == nextpid)
   root = true, foreach_child_break(); end_foreach_child(); }
      if (root && cell.pid != nextpid) {
  { foreach_neighbor()
   if (cell.pid != nextpid && is_newpid()) {
     if (fp)
       fprintf (fp, "%g %g %g %d %d root\n",
         x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
     cell.flags |= sent;
   } end_foreach_neighbor(); }
      }
    }

    if ((is_local(cell) && ((NewPid *)&val(newpid,0,0,0))->pid - 1 == nextpid) || root) {
       { foreach_neighbor(1)
 if (cell.neighbors && cell.pid != nextpid)
    { foreach_child()
     if (cell.pid != nextpid && is_newpid()) {
       if (fp)
  fprintf (fp, "%g %g %g %d %d nextpid\n",
    x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
       cell.flags |= sent;
     } end_foreach_child(); } end_foreach_neighbor(); }
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  return tree (sizeof(Cell) + datasize, newpid);
}

static void send_tree (Array * a, int to, MPI_Request * r)
{
  MPI_Isend (&a->len, 1, MPI_LONG, to, (256), MPI_COMM_WORLD, &r[0]);
  if (a->len > 0) {
    MPI_Isend (a->p, a->len, MPI_BYTE, to, (256), MPI_COMM_WORLD, &r[1]);
    ((Quadtree *)grid)->dirty = true;
  }
}

static void receive_tree (int from, scalar newpid, FILE * fp)
{
  Array a;
  mpi_recv_check (&a.len, 1, MPI_LONG, from, (256),
    MPI_COMM_WORLD, MPI_STATUS_IGNORE, "receive_tree (len)");
  if (a.len > 0) {
    a.p = pmalloc (a.len,__func__,__FILE__,__LINE__);
    if (fp)
      fprintf (fp, "receiving %ld from %d\n", a.len, from);
    mpi_recv_check (a.p, a.len, MPI_BYTE, from, (256),
      MPI_COMM_WORLD, MPI_STATUS_IGNORE, "receive_tree (p)");

     { foreach_tree (&a, sizeof(Cell) + datasize, NULL){

#line 156 "/home/popinet/basilisk-octree/src/grid/balance.h"
 {
      memcpy (((char *)&cell) + sizeof(Cell), ((char *)c) + sizeof(Cell),
       datasize);
      assert (((NewPid *)&val(newpid,0,0,0))->pid > 0);
      if (fp)
 fprintf (fp, "%g %g %g %d %d %d %d %d %d recv\n",
   x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid,
   c->flags & leaf,
   cell.flags & leaf, from, ((NewPid *)&val(newpid,0,0,0))->leaf);
    } } end_foreach_tree(); }
    pfree (a.p,__func__,__FILE__,__LINE__);
    ((Quadtree *)grid)->dirty = true;
  }
}

static void wait_tree (Array * a, MPI_Request * r)
{
  MPI_Wait (&r[0], MPI_STATUS_IGNORE);
  if (a->len > 0)
    MPI_Wait (&r[1], MPI_STATUS_IGNORE);
}

static void check_flags()
{







}

struct {
  int min;
  bool leaves;

  int npe;
} mpi = {
  10000,
  true
};


bool balance()
{ trace ("balance", "/home/popinet/basilisk-octree/src/grid/balance.h", 201);
  if (npe() == 1)
    { bool _ret =  false; end_trace("balance", "/home/popinet/basilisk-octree/src/grid/balance.h", 203);  return _ret; }

  assert (sizeof(NewPid) == sizeof(double));

  check_flags();

  long nl = 0;

  if (mpi.leaves)
     { foreach(){

#line 212 "/home/popinet/basilisk-octree/src/grid/balance.h"

      nl++; } end_foreach(); }
  else
     { foreach_cell(){

#line 215 "/home/popinet/basilisk-octree/src/grid/balance.h"
 {
      if (is_local(cell))
 nl++;
      if (is_leaf(cell))
 continue;
    } } end_foreach_cell(); }

  long nt = nl, nmin = nl, nmax = nl;

  mpi_all_reduce (nt, MPI_LONG, MPI_SUM);
  mpi_all_reduce (nmax, MPI_LONG, MPI_MAX);
  mpi_all_reduce (nmin, MPI_LONG, MPI_MIN);
  long ne = max(1, nt/npe());

  if (ne < mpi.min) {
    mpi.npe = max(1, nt/mpi.min);
    ne = max(1, nt/mpi.npe);
  }
  else
    mpi.npe = npe();

  if (nmax - nmin <= 1)
    { bool _ret =  false; end_trace("balance", "/home/popinet/basilisk-octree/src/grid/balance.h", 237);  return _ret; }

  scalar newpid= new_scalar("newpid");
  double zn = z_indexing (newpid, mpi.leaves);
  if (pid() == 0)
    assert (zn + 1 == nt);

  FILE * fp = NULL;
#line 254 "/home/popinet/basilisk-octree/src/grid/balance.h"
  bool next = false, prev = false;
   { foreach_cell_all(){

#line 255 "/home/popinet/basilisk-octree/src/grid/balance.h"
 {
    if (is_local(cell)) {
      int pid = balanced_pid (val(newpid,0,0,0), nt, mpi.npe);
      pid = clamp (pid, cell.pid - 1, cell.pid + 1);
      if (pid == pid() + 1)
 next = true;
      else if (pid == pid() - 1)
 prev = true;
      ((NewPid *)&val(newpid,0,0,0))->pid = pid + 1;
      ((NewPid *)&val(newpid,0,0,0))->leaf = is_leaf(cell);
      ((NewPid *)&val(newpid,0,0,0))->prolongation = (!is_leaf(cell) && !cell.neighbors && cell.pid >= 0);
      if (fp)
 fprintf (fp, "%g %g %d %d newpid\n", x, y, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
    }
    else
      val(newpid,0,0,0) = 0;
  } } end_foreach_cell_all(); }
  for (int l = 0; l <= depth(); l++)
    { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){newpid,{-1}}), l); };
#line 297 "/home/popinet/basilisk-octree/src/grid/balance.h"
  Array * anext = next ? neighborhood (newpid, pid() + 1, fp) : array_new();
  Array * aprev = prev ? neighborhood (newpid, pid() - 1, fp) : array_new();

  if (fp)
    fflush (fp);

  check_flags();


  MPI_Request rprev[2], rnext[2];
  if (pid() > 0)
    send_tree (aprev, pid() - 1, rprev);
  if (pid() < npe() - 1)
    send_tree (anext, pid() + 1, rnext);


  if (pid() < npe() - 1)
    receive_tree (pid() + 1, newpid, fp);
  if (pid() > 0)
    receive_tree (pid() - 1, newpid, fp);


  if (pid() > 0)
    wait_tree (aprev, rprev);
  array_free (aprev);
  if (pid() < npe() - 1)
    wait_tree (anext, rnext);
  array_free (anext);

  if (fp)
    fflush (fp);


  int pid_changed = false;
   { foreach_cell_all(){

#line 331 "/home/popinet/basilisk-octree/src/grid/balance.h"
 {
    if (cell.pid >= 0) {
      if (is_newpid()) {
 if (fp)
   fprintf (fp, "%g %g %g %d %d %d %d %d new\n",
     x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid,
     is_leaf(cell), cell.neighbors, ((NewPid *)&val(newpid,0,0,0))->leaf);
 if (cell.pid != ((NewPid *)&val(newpid,0,0,0))->pid - 1) {
   cell.pid = ((NewPid *)&val(newpid,0,0,0))->pid - 1;
   cell.flags &= ~(active|border);
   if (is_local(cell))
     cell.flags |= active;
   pid_changed = true;
 }
 if (((NewPid *)&val(newpid,0,0,0))->leaf && !is_leaf(cell) && cell.neighbors)
   coarsen_cell_recursive (point, NULL);
      }
      else if (level > 0 && ((NewPid *)&coarse(newpid,0,0,0))->leaf)
 cell.pid = aparent(0,0,0).pid;
    }

    if (!cell.neighbors && allocated_child(0,0,0)) {
      if (fp)
 fprintf (fp, "%g %g %g %d %d freechildren\n",
   x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
      free_children (point);
    }
  } } end_foreach_cell_all(); }

  if (((Quadtree *)grid)->dirty || pid_changed) {


     { foreach_cell_post (!is_leaf (cell)){

#line 363 "/home/popinet/basilisk-octree/src/grid/balance.h"

      if (!is_leaf(cell) && !is_local(cell)) {
 unsigned short flags = cell.flags & ~active;
  { foreach_child()
   if (is_active(cell))
     flags |= active, foreach_child_break(); end_foreach_child(); }
 cell.flags = flags;
      } } end_foreach_cell_post(); }

    flag_border_cells();
    pid_changed = true;
  }

  if (fp)
    fclose (fp);

  mpi_all_reduce (pid_changed, MPI_INT, MPI_MAX);
  if (pid_changed)
    mpi_boundary_update();

  { bool _ret =  pid_changed; delete (((scalar []){newpid,{-1}}));  end_trace("balance", "/home/popinet/basilisk-octree/src/grid/balance.h", 383);  return _ret; }
 delete (((scalar []){newpid,{-1}}));  end_trace("balance", "/home/popinet/basilisk-octree/src/grid/balance.h", 384); }
#line 1717 "/home/popinet/basilisk-octree/src/grid/tree.h"
#endif
#line 4 "/home/popinet/basilisk-octree/src/grid/octree.h"

void octree_methods() {
  quadtree_methods();
}
#line 15 "atomisation-cpp.c"
static double _boundary0 (Point point, Point neighbor, scalar _s);
static double _boundary0_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary1 (Point point, Point neighbor, scalar _s);
static double _boundary1_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary2 (Point point, Point neighbor, scalar _s);
static double _boundary2_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary3 (Point point, Point neighbor, scalar _s);
static double _boundary3_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary4 (Point point, Point neighbor, scalar _s);
static double _boundary4_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary5 (Point point, Point neighbor, scalar _s);
static double _boundary5_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary6 (Point point, Point neighbor, scalar _s);
static double _boundary6_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary7 (Point point, Point neighbor, scalar _s);
static double _boundary7_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary8 (Point point, Point neighbor, scalar _s);
static double _boundary8_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary9 (Point point, Point neighbor, scalar _s);
static double _boundary9_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary10 (Point point, Point neighbor, scalar _s);
static double _boundary10_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary11 (Point point, Point neighbor, scalar _s);
static double _boundary11_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary12 (Point point, Point neighbor, scalar _s);
static double _boundary12_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary13 (Point point, Point neighbor, scalar _s);
static double _boundary13_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary14 (Point point, Point neighbor, scalar _s);
static double _boundary14_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary15 (Point point, Point neighbor, scalar _s);
static double _boundary15_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary16 (Point point, Point neighbor, scalar _s);
static double _boundary16_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary17 (Point point, Point neighbor, scalar _s);
static double _boundary17_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary18 (Point point, Point neighbor, scalar _s);
static double _boundary18_homogeneous (Point point, Point neighbor, scalar _s);
static double _boundary19 (Point point, Point neighbor, scalar _s);
static double _boundary19_homogeneous (Point point, Point neighbor, scalar _s);
#line 1 "atomisation.c"

#line 1 "navier-stokes/centered.h"
#line 1 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
#line 24 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
#line 1 "./run.h"
#line 1 "/home/popinet/basilisk-octree/src/run.h"






#line 1 "./utils.h"
#line 1 "/home/popinet/basilisk-octree/src/utils.h"


double DT = 1e10;

double CFL = 0.5;

void runge_kutta (int stages,
    double t, double dt,
    int nv, scalar f[nv], scalar df[stages][nv],
    void (* advance) (double t, scalar f[nv], scalar df[nv]),
    void (* update) (double t, scalar f[nv]))
{
  switch (stages) {
  case 1:
    (* advance) (t, f, df[0]);
     { foreach(){

#line 16 "/home/popinet/basilisk-octree/src/utils.h"

      for (int v = 0; v < nv; v++)
 val(f[v],0,0,0) += val(df[0][v],0,0,0)*dt; } end_foreach(); }
    (* update) (t + dt, f);
    break;

  case 2:
    (* advance) (t, f, df[0]);
     { foreach(){

#line 24 "/home/popinet/basilisk-octree/src/utils.h"

      for (int v = 0; v < nv; v++)
 val(df[0][v],0,0,0) = val(f[v],0,0,0) + val(df[0][v],0,0,0)*dt/2.; } end_foreach(); }
    (* update) (t + dt/2., df[0]);

    (* advance) (t + dt/2., df[0], df[1]);
     { foreach(){

#line 30 "/home/popinet/basilisk-octree/src/utils.h"

      for (int v = 0; v < nv; v++)
 val(f[v],0,0,0) += val(df[1][v],0,0,0)*dt; } end_foreach(); }
    (* update) (t + dt, f);
    break;

  default:

    assert(false);
  }
}

double change (scalar v, scalar vn)
{
  double max = 0.;
   { 
#undef _OMPSTART
#define _OMPSTART double _max = max; 
#undef _OMPEND
#define _OMPEND OMP(omp critical) if (_max > max) max = _max; mpi_all_reduce_double (max, MPI_MAX); 
#line 45
foreach(){

#line 45 "/home/popinet/basilisk-octree/src/utils.h"
 {
    double dv = fabs (val(v,0,0,0) - val(vn,0,0,0));
    if (dv > _max)
      _max = dv;
    val(vn,0,0,0) = val(v,0,0,0);
  } } end_foreach();
#undef _OMPSTART
#undef _OMPEND
#define _OMPSTART
#define _OMPEND
#line 51
 }
  return max;
}

typedef struct {
  double cpu;
  double real;
  double speed;
  double min;
  double avg;
  double max;
  size_t tnc;
  long mem;
} timing;

timing timer_timing (timer t, int i, size_t tnc, double * mpi)
{
  timing s;
#if _MPI
  s.avg = mpi_time - t.tm;
#endif
  clock_t end = clock();
  s.cpu = ((double) (end - t.c))/CLOCKS_PER_SEC;
  s.real = timer_elapsed (t);
  if (tnc == 0) {
     { 
#undef _OMPSTART
#define _OMPSTART 
#undef _OMPEND
#define _OMPEND mpi_all_reduce_double (tnc, MPI_SUM); 
#line 75
foreach(reduction(+:tnc)){

#line 75 "/home/popinet/basilisk-octree/src/utils.h"
 tnc++; } end_foreach();
#undef _OMPSTART
#undef _OMPEND
#define _OMPSTART
#define _OMPEND
#line 76
 }
    s.tnc = tnc;
    tnc *= i;
  }
  else
    s.tnc = tnc;
#if (_GNU_SOURCE || _DARWIN_C_SOURCE)
  struct rusage usage;
  getrusage (RUSAGE_SELF, &usage);
  s.mem = usage.ru_maxrss;
#else
  s.mem = 0;
#endif
#if _MPI
  if (mpi)
    MPI_Allgather (&s.avg, 1, MPI_DOUBLE, mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
  s.max = s.min = s.avg;
  mpi_all_reduce (s.max, MPI_DOUBLE, MPI_MAX);
  mpi_all_reduce (s.min, MPI_DOUBLE, MPI_MIN);
  mpi_all_reduce (s.avg, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.real, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.mem, MPI_LONG, MPI_SUM);
  s.real /= npe();
  s.avg /= npe();
  s.mem /= npe();
#else
  s.min = s.max = s.avg = 0.;
#endif
  s.speed = s.real > 0. ? tnc/s.real : -1;
  return s;
}

void timer_print (timer t, int i, size_t tnc)
{
  timing s = timer_timing (t, i, tnc, NULL);
  fprintf (fout,
    "\n# " "Octree"
    ", %d steps, %g CPU, %.4g real, %.3g points.step/s, %d var\n",
    i, s.cpu, s.real, s.speed, (int) (datasize/sizeof(double)));
#if _MPI
  fprintf (fout,
    "# %d procs, MPI: min %.2g (%.2g%%) "
    "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
    npe(),
    s.min, 100.*s.min/s.real,
    s.avg, 100.*s.avg/s.real,
    s.max, 100.*s.max/s.real);
#endif
}

typedef struct {
  double avg, rms, max, volume;
} norm;

norm normf (scalar f)
{
  double avg = 0., rms = 0., max = 0., volume = 0.;
   { 
#undef _OMPSTART
#define _OMPSTART double _max = max; 
#undef _OMPEND
#define _OMPEND OMP(omp critical) if (_max > max) max = _max; mpi_all_reduce_double (max, MPI_MAX); mpi_all_reduce_double (avg, MPI_SUM); mpi_all_reduce_double (rms, MPI_SUM); mpi_all_reduce_double (volume, MPI_SUM); 
#line 132

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 132
foreach( reduction(+:avg)
   reduction(+:rms) reduction(+:volume)){

#line 133 "/home/popinet/basilisk-octree/src/utils.h"

    if (val(f,0,0,0) != nodata) {
      double v = fabs(val(f,0,0,0));
      if (v > _max) _max = v;
      volume += (cube(Delta)*val_cm(cm,0,0,0));
      avg += (cube(Delta)*val_cm(cm,0,0,0))*v;
      rms += (cube(Delta)*val_cm(cm,0,0,0))*sq(v);
    } } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 132
foreach( reduction(+:avg)
   reduction(+:rms) reduction(+:volume)){

#line 133 "/home/popinet/basilisk-octree/src/utils.h"

    if (val(f,0,0,0) != nodata) {
      double v = fabs(val(f,0,0,0));
      if (v > _max) _max = v;
      volume += (cube(Delta)*val_cm(cm,0,0,0));
      avg += (cube(Delta)*val_cm(cm,0,0,0))*v;
      rms += (cube(Delta)*val_cm(cm,0,0,0))*sq(v);
    } } end_foreach(); }
#undef _OMPSTART
#undef _OMPEND
#define _OMPSTART
#define _OMPEND
#line 141
 }
  norm n;
  n.avg = avg/volume;
  n.rms = sqrt(rms/volume);
  n.max = max;
  n.volume = volume;
  return n;
}

typedef struct {
  double min, max, sum, stddev, volume;
} stats;

stats statsf (scalar f)
{
  double min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
   { 
#undef _OMPSTART
#define _OMPSTART double _max = max; double _min = min; 
#undef _OMPEND
#define _OMPEND mpi_all_reduce_double (sum, MPI_SUM); mpi_all_reduce_double (sum2, MPI_SUM); mpi_all_reduce_double (volume, MPI_SUM); OMP(omp critical) if (_max > max) max = _max; mpi_all_reduce_double (max, MPI_MAX); OMP(omp critical) if (_min < min) min = _min; mpi_all_reduce_double (min, MPI_MIN); 
#line 156

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 156
foreach(reduction(+:sum) reduction(+:sum2) reduction(+:volume)){

#line 157 "/home/popinet/basilisk-octree/src/utils.h"

    if (val(f,0,0,0) != nodata) {
      volume += (cube(Delta)*val_cm(cm,0,0,0));
      sum += (cube(Delta)*val_cm(cm,0,0,0))*val(f,0,0,0);
      sum2 += (cube(Delta)*val_cm(cm,0,0,0))*sq(val(f,0,0,0));
      if (val(f,0,0,0) > _max) _max = val(f,0,0,0);
      if (val(f,0,0,0) < _min) _min = val(f,0,0,0);
    } } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 156
foreach(reduction(+:sum) reduction(+:sum2) reduction(+:volume)){

#line 157 "/home/popinet/basilisk-octree/src/utils.h"

    if (val(f,0,0,0) != nodata) {
      volume += (cube(Delta)*val_cm(cm,0,0,0));
      sum += (cube(Delta)*val_cm(cm,0,0,0))*val(f,0,0,0);
      sum2 += (cube(Delta)*val_cm(cm,0,0,0))*sq(val(f,0,0,0));
      if (val(f,0,0,0) > _max) _max = val(f,0,0,0);
      if (val(f,0,0,0) < _min) _min = val(f,0,0,0);
    } } end_foreach(); }
#undef _OMPSTART
#undef _OMPEND
#define _OMPSTART
#define _OMPEND
#line 165
 }
  stats s;
  s.min = min, s.max = max, s.sum = sum, s.volume = volume;
  sum2 -= sum*sum/volume;
  s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
  return s;
}

static double generic_limiter (double r, double beta)
{
  double v1 = min (r, beta), v2 = min (beta*r, 1.);
  v1 = max (0., v1);
  return max (v1, v2);
}

double minmod (double s0, double s1, double s2) {
  return generic_limiter ((s2 - s1)/(s1 - s0), 1.)*(s1 - s0);
}

double superbee (double s0, double s1, double s2) {
  return generic_limiter ((s2 - s1)/(s1 - s0), 2.)*(s1 - s0);
}

double sweby (double s0, double s1, double s2) {
  return generic_limiter ((s2 - s1)/(s1 - s0), 1.5)*(s1 - s0);
}

double theta = 1.3;

double minmod2 (double s0, double s1, double s2)
{
  if (s0 < s1 && s1 < s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 < d1) d1 = d2;
    return min(d1, d3);
  }
  if (s0 > s1 && s1 > s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 > d1) d1 = d2;
    return max(d1, d3);
  }
  return 0.;
}

void gradients (scalar * f, vector * g)
{
  assert (list_len(f) == vectors_len(g));
   { foreach(){

#line 211 "/home/popinet/basilisk-octree/src/utils.h"
 {
    scalar s; vector v;
    scalar * _i0 = f; vector * _i1 = g; if (f) for (s = *f, v = *g; ((scalar *)&s)->i >= 0; s = *++_i0, v = *++_i1) {
      if (_attribute[s.i].gradient)
 {
#line 215

   val(v.x,0,0,0) = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0))/Delta;
#line 215

   val(v.y,0,0,0) = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0))/Delta;
#line 215

   val(v.z,0,0,0) = _attribute[s.i].gradient (val(s,0,0,-1), val(s,0,0,0), val(s,0,0,1))/Delta;}
      else
 {
#line 218

   val(v.x,0,0,0) = (val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta);
#line 218

   val(v.y,0,0,0) = (val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta);
#line 218

   val(v.z,0,0,0) = (val(s,0,0,1) - val(s,0,0,-1))/(2.*Delta);}
    }
  } } end_foreach(); }
  boundary ((scalar *) g);
}
#line 235 "/home/popinet/basilisk-octree/src/utils.h"
void vorticity (const vector u, scalar omega)
{

     { foreach(){

#line 238 "/home/popinet/basilisk-octree/src/utils.h"

      val(omega,0,0,0) = (val(u.y,1,0,0) - val(u.y,-1,0,0) + val(u.x,0,-1,0) - val(u.x,0,1,0))/(2.*Delta); } end_foreach(); }
    boundary (((scalar []){omega,{-1}}));

}




struct {

  long n;

  long tn;

  long nc;

  long tnc;

  double t;

  double speed;

  timer gt;
} perf;

static void update_perf() {
  long tn = 0;
  perf.n = 0;
   { 
#undef _OMPSTART
#define _OMPSTART 
#undef _OMPEND
#define _OMPEND mpi_all_reduce_double (tn, MPI_SUM); 
#line 267
foreach(reduction(+:tn)){

#line 267 "/home/popinet/basilisk-octree/src/utils.h"

    perf.n++, tn++; } end_foreach();
#undef _OMPSTART
#undef _OMPEND
#define _OMPSTART
#define _OMPEND
#line 269
 }
  perf.tn = tn;
  perf.nc += perf.n;
  perf.tnc += perf.tn;
  perf.t = timer_elapsed (perf.gt);
  perf.speed = perf.tnc/perf.t;
}

#line 1 "./output.h"
#line 1 "/home/popinet/basilisk-octree/src/output.h"
#line 33 "/home/popinet/basilisk-octree/src/output.h"
struct OutputField {
  scalar * list;
  FILE * fp;
  int n;
  bool linear;
};


void output_field (struct OutputField p)
{ trace ("output_field", "/home/popinet/basilisk-octree/src/output.h", 42);
  if (!p.list) p.list = all;
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = qstdout();

  int len = list_len(p.list);
  double ** field = matrix_new (p.n, p.n, len*sizeof(double));

  double Delta = L0/p.n;
  for (int i = 0; i < p.n; i++) {
    double xp = Delta*i + X0 + Delta/2.;
    for (int j = 0; j < p.n; j++) {
      double yp = Delta*j + Y0 + Delta/2.;
      if (p.linear) {
 int k = 0;
 if (p.list) for (scalar s = *p.list, *_i63 = p.list; ((scalar *)&s)->i >= 0; s = *++_i63)
   field[i][len*j + k++] = interpolate ((struct _interpolate){s, xp, yp});
      }
      else {
 Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 61 "/home/popinet/basilisk-octree/src/output.h"

 int k = 0;
 if (p.list) for (scalar s = *p.list, *_i64 = p.list; ((scalar *)&s)->i >= 0; s = *++_i64)
   field[i][len*j + k++] = point.level >= 0 ? val(s,0,0,0) : nodata;
      }
    }
  }

  if (pid() == 0) {
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], len*p.n*p.n, MPI_DOUBLE, MPI_MIN, 0,
  MPI_COMM_WORLD);
#endif
    fprintf (p.fp, "# 1:x 2:y");
    int i = 3;
    if (p.list) for (scalar s = *p.list, *_i65 = p.list; ((scalar *)&s)->i >= 0; s = *++_i65)
      fprintf (p.fp, " %d:%s", i++, _attribute[s.i].name);
    fputc('\n', p.fp);
    for (int i = 0; i < p.n; i++) {
      double xp = Delta*i + X0 + Delta/2.;
      for (int j = 0; j < p.n; j++) {
 double yp = Delta*j + Y0 + Delta/2.;
 fprintf (p.fp, "%g %g", xp, yp);
 int k = 0;
 if (p.list) for (scalar s = *p.list, *_i66 = p.list; ((scalar *)&s)->i >= 0; s = *++_i66)
   fprintf (p.fp, " %g", field[i][len*j + k++]);
 fputc ('\n', p.fp);
      }
      fputc ('\n', p.fp);
    }
    fflush (p.fp);
  }
#if _MPI
  else
    MPI_Reduce (field[0], NULL, len*p.n*p.n, MPI_DOUBLE, MPI_MIN, 0,
  MPI_COMM_WORLD);
#endif

  matrix_free (field);
 end_trace("output_field", "/home/popinet/basilisk-octree/src/output.h", 100); }
#line 128 "/home/popinet/basilisk-octree/src/output.h"
struct OutputMatrix {
  scalar f;
  FILE * fp;
  int n;
  bool linear;
};


void output_matrix (struct OutputMatrix p)
{ trace ("output_matrix", "/home/popinet/basilisk-octree/src/output.h", 137);
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = qstdout();
  float fn = p.n;
  float Delta = L0/fn;
  fwrite (&fn, sizeof(float), 1, p.fp);
  for (int j = 0; j < p.n; j++) {
    float yp = Delta*j + X0 + Delta/2.;
    fwrite (&yp, sizeof(float), 1, p.fp);
  }
  for (int i = 0; i < p.n; i++) {
    float xp = Delta*i + X0 + Delta/2.;
    fwrite (&xp, sizeof(float), 1, p.fp);
    for (int j = 0; j < p.n; j++) {
      float yp = Delta*j + Y0 + Delta/2., v;
      if (p.linear)
 v = interpolate ((struct _interpolate){p.f, xp, yp});
      else {
 Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 155 "/home/popinet/basilisk-octree/src/output.h"

 assert (point.level >= 0);
 v = val(p.f,0,0,0);
      }
      fwrite (&v, sizeof(float), 1, p.fp);
    }
  }
  fflush (p.fp);
 end_trace("output_matrix", "/home/popinet/basilisk-octree/src/output.h", 163); }
#line 172 "/home/popinet/basilisk-octree/src/output.h"
typedef void (* colormap) (double cmap[127][3]);

void jet (double cmap[127][3])
{
  for (int i = 0; i < 127; i++) {
    cmap[i][0] =
      i <= 46 ? 0. :
      i >= 111 ? -0.03125*(i - 111) + 1. :
      i >= 78 ? 1. :
      0.03125*(i - 46);
    cmap[i][1] =
      i <= 14 || i >= 111 ? 0. :
      i >= 79 ? -0.03125*(i - 111) :
      i <= 46 ? 0.03125*(i - 14) :
      1.;
    cmap[i][2] =
      i >= 79 ? 0. :
      i >= 47 ? -0.03125*(i - 79) :
      i <= 14 ? 0.03125*(i - 14) + 1.:
      1.;
  }
}

void cool_warm (double cmap[127][3])
{






  static double basemap[33][3] = {
    {0.2298057, 0.298717966, 0.753683153},
    {0.26623388, 0.353094838, 0.801466763},
    {0.30386891, 0.406535296, 0.84495867},
    {0.342804478, 0.458757618, 0.883725899},
    {0.38301334, 0.50941904, 0.917387822},
    {0.424369608, 0.558148092, 0.945619588},
    {0.46666708, 0.604562568, 0.968154911},
    {0.509635204, 0.648280772, 0.98478814},
    {0.552953156, 0.688929332, 0.995375608},
    {0.596262162, 0.726149107, 0.999836203},
    {0.639176211, 0.759599947, 0.998151185},
    {0.681291281, 0.788964712, 0.990363227},
    {0.722193294, 0.813952739, 0.976574709},
    {0.761464949, 0.834302879, 0.956945269},
    {0.798691636, 0.849786142, 0.931688648},
    {0.833466556, 0.860207984, 0.901068838},
    {0.865395197, 0.86541021, 0.865395561},
    {0.897787179, 0.848937047, 0.820880546},
    {0.924127593, 0.827384882, 0.774508472},
    {0.944468518, 0.800927443, 0.726736146},
    {0.958852946, 0.769767752, 0.678007945},
    {0.96732803, 0.734132809, 0.628751763},
    {0.969954137, 0.694266682, 0.579375448},
    {0.966811177, 0.650421156, 0.530263762},
    {0.958003065, 0.602842431, 0.481775914},
    {0.943660866, 0.551750968, 0.434243684},
    {0.923944917, 0.49730856, 0.387970225},
    {0.89904617, 0.439559467, 0.343229596},
    {0.869186849, 0.378313092, 0.300267182},
    {0.834620542, 0.312874446, 0.259301199},
    {0.795631745, 0.24128379, 0.220525627},
    {0.752534934, 0.157246067, 0.184115123},
    {0.705673158, 0.01555616, 0.150232812}
  };

  for (int i = 0; i < 127; i++) {
    double x = i*(32 - 1e-10)/(127 - 1);
    int j = x; x -= j;
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (1. - x)*basemap[j][k] + x*basemap[j+1][k];
  }
}

void gray (double cmap[127][3])
{
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = i/(127 - 1.);
}

void randomap (double cmap[127][3])
{
  srand(0);
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (noise() + 1.)/2.;
}





typedef struct {
  unsigned char r, g, b;
} color;

color colormap_color (double cmap[127][3],
        double val, double min, double max)
{
  color c;
  if (val == nodata) {
    c.r = c.g = c.b = 0;
    return c;
  }
  val = val <= min ? 0. : val >= max ? 0.9999 : (val - min)/(max - min);
  int i = val*(127 - 1);
  double coef = val*(127 - 1) - i;
  assert (i < 127 - 1);
  unsigned char * c1 = (unsigned char *) &c;
  for (int j = 0; j < 3; j++)
    c1[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);
  return c;
}
#line 347 "/home/popinet/basilisk-octree/src/output.h"
struct OutputPPM {
  scalar f;
  FILE * fp;
  int n;
  char * file;
  double min, max, spread, z;
  bool linear;
  double box[2][2];
  scalar mask;
  colormap map;
};


void output_ppm (struct OutputPPM p)
{ trace ("output_ppm", "/home/popinet/basilisk-octree/src/output.h", 361);

  if (p.n == 0) p.n = N;
  if (p.min == 0 && p.max == 0) {
    stats s = statsf (p.f);
    double avg = s.sum/s.volume, spread = (p.spread ? p.spread : 5.)*s.stddev;
    p.min = avg - spread; p.max = avg + spread;
  }
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }
  if (!p.map)
    p.map = jet;

  double fn = p.n;
  double Delta = (p.box[1][0] - p.box[0][0])/fn;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;

  color ** ppm = matrix_new (ny, p.n, sizeof(color));
  double cmap[127][3];
  p.map (cmap);
  OMP_PARALLEL()
  OMP(omp for schedule(static))
  for (int j = 0; j < ny; j++) {
    double yp = Delta*j + p.box[0][1] + Delta/2.;
    for (int i = 0; i < p.n; i++) {
      double xp = Delta*i + p.box[0][0] + Delta/2., v;
      if (p.mask.i) {
 if (p.linear) {
   double m = interpolate ((struct _interpolate){p.mask, xp, yp, p.z});
   if (m < 0.)
     v = nodata;
   else
     v = interpolate ((struct _interpolate){p.f, xp, yp, p.z});
 }
 else {
   Point point = locate ((struct _locate){xp, yp, p.z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 399 "/home/popinet/basilisk-octree/src/output.h"

   if (point.level < 0 || val(p.mask,0,0,0) < 0.)
     v = nodata;
   else
     v = val(p.f,0,0,0);
 }
      }
      else if (p.linear)
 v = interpolate ((struct _interpolate){p.f, xp, yp, p.z});
      else {
 Point point = locate ((struct _locate){xp, yp, p.z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 409 "/home/popinet/basilisk-octree/src/output.h"

 v = point.level >= 0 ? val(p.f,0,0,0) : nodata;
      }
      ppm[ny - 1 - j][i] = colormap_color (cmap, v, p.min, p.max);
    }
  }
  OMP_END_PARALLEL()

  if (pid() == 0) {
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, ppm[0], 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
  MPI_COMM_WORLD);
#endif
    if (!p.fp) p.fp = qstdout();
    if (p.file) {
      char * command = pmalloc (strlen ("convert ppm:- ") + strlen (p.file) + 1,__func__,__FILE__,__LINE__);
      strcpy (command, "convert ppm:- ");
      strcat (command, p.file);
      p.fp = qpopen (command, "w");
      pfree (command,__func__,__FILE__,__LINE__);
    }

    fprintf (p.fp, "P6\n%u %u 255\n", p.n, ny);
    fwrite (((void **) ppm)[0], sizeof(color), ny*p.n, p.fp);

    if (p.file)
      qpclose (p.fp);
    else
      fflush (p.fp);
  }
#if _MPI
  else
    MPI_Reduce (ppm[0], NULL, 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
  MPI_COMM_WORLD);
#endif

  matrix_free (ppm);
 end_trace("output_ppm", "/home/popinet/basilisk-octree/src/output.h", 446); }
#line 487 "/home/popinet/basilisk-octree/src/output.h"
struct OutputGRD {
  scalar f;
  FILE * fp;
  double Delta;
  bool linear;
  double box[2][2];
  scalar mask;
};


void output_grd (struct OutputGRD p)
{ trace ("output_grd", "/home/popinet/basilisk-octree/src/output.h", 498);

  if (!p.fp) p.fp = qstdout();
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
    if (p.Delta == 0) p.Delta = L0/N;
  }

  double Delta = p.Delta;
  int nx = (p.box[1][0] - p.box[0][0])/Delta;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;


  fprintf (p.fp, "ncols          %d\n", nx);
  fprintf (p.fp, "nrows          %d\n", ny);
  fprintf (p.fp, "xllcorner      %g\n", p.box[0][0]);
  fprintf (p.fp, "yllcorner      %g\n", p.box[0][1]);
  fprintf (p.fp, "cellsize       %g\n", Delta);
  fprintf (p.fp, "nodata_value   -9999\n");


  for (int j = ny-1; j >= 0; j--) {
    double yp = Delta*j + p.box[0][1] + Delta/2.;
    for (int i = 0; i < nx; i++) {
      double xp = Delta*i + p.box[0][0] + Delta/2., v;
      if (p.mask.i) {
 if (p.linear) {
   double m = interpolate ((struct _interpolate){p.mask, xp, yp});
   if (m < 0.)
     v = nodata;
   else
     v = interpolate ((struct _interpolate){p.f, xp, yp});
 }
 else {
   Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 534 "/home/popinet/basilisk-octree/src/output.h"

   if (point.level < 0 || val(p.mask,0,0,0) < 0.)
     v = nodata;
   else
     v = val(p.f,0,0,0);
 }
      }
      else if (p.linear)
 v = interpolate ((struct _interpolate){p.f, xp, yp});
      else {
 Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 544 "/home/popinet/basilisk-octree/src/output.h"

 v = point.level >= 0 ? val(p.f,0,0,0) : nodata;
      }
      if (v == nodata)
 fprintf (p.fp, "-9999 ");
      else
 fprintf (p.fp, "%f ", v);
    }
    fprintf (p.fp, "\n");
  }

  fflush (p.fp);
 end_trace("output_grd", "/home/popinet/basilisk-octree/src/output.h", 556); }
#line 586 "/home/popinet/basilisk-octree/src/output.h"
struct OutputGfs {
  FILE * fp;
  scalar * list;
  double t;
  char * file;
  bool translate;
};

static char * replace (const char * input, int target, int with,
         bool translate)
{
  if (translate) {
    if (!strcmp (input, "u.x"))
      return pstrdup ("U",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.y"))
      return pstrdup ("V",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.z"))
      return pstrdup ("W",__func__,__FILE__,__LINE__);
  }
  char * name = pstrdup (input,__func__,__FILE__,__LINE__), * i = name;
  while (*i != '\0') {
    if (*i == target)
      *i = with;
    i++;
  }
  return name;
}


void output_gfs (struct OutputGfs p)
{ trace ("output_gfs", "/home/popinet/basilisk-octree/src/output.h", 616);
  bool opened = false;
  if (p.fp == NULL) {
    if (p.file == NULL)
      p.fp = qstdout();
    else if (!(p.fp = fopen (p.file, "w"))) {
      perror (p.file);
      exit (1);
    }
    else
      opened = true;
  }
  scalar * list = p.list ? p.list : list_copy (all);

  fprintf (p.fp,
    "1 0 GfsSimulation GfsBox GfsGEdge { binary = 1"
    " x = %g y = %g ",
    0.5 + X0/L0, 0.5 + Y0/L0);

  fprintf (p.fp, "z = %g ", 0.5 + Z0/L0);


  if (list != NULL && list[0].i != -1) {
    scalar s = list[0];
    char * name = replace (_attribute[s.i].name, '.', '_', p.translate);
    fprintf (p.fp, "variables = %s", name);
    pfree (name,__func__,__FILE__,__LINE__);
    for (int i = 1; i < list_len(list); i++) {
      scalar s = list[i];
      if (_attribute[s.i].name) {
 char * name = replace (_attribute[s.i].name, '.', '_', p.translate);
 fprintf (p.fp, ",%s", name);
 pfree (name,__func__,__FILE__,__LINE__);
      }
    }
    fprintf (p.fp, " ");
  }
  fprintf (p.fp, "} {\n");
  if (p.t > 0.)
    fprintf (p.fp, "  Time { t = %g }\n", p.t);
  if (L0 != 1.)
    fprintf (p.fp, "  PhysicalParams { L = %g }\n", L0);
  fputs ("  VariableTracerVOF f\n", p.fp);
  fprintf (p.fp, "}\nGfsBox { x = 0 y = 0 z = 0 } {\n");

#if _MPI
  long header;
  if ((header = ftell (p.fp)) < 0) {
    perror ("output_gfs(): error in header");
    exit (1);
  }
  int cell_size = sizeof(unsigned) + sizeof(double);
  if (list) for (scalar s = *list, *_i67 = list; ((scalar *)&s)->i >= 0; s = *++_i67)
    if (_attribute[s.i].name)
      cell_size += sizeof(double);
  scalar index = new_scalar("index");
  size_t total_size = header + (z_indexing (index, false) + 1)*cell_size;
#endif



   { foreach_cell(){

#line 677 "/home/popinet/basilisk-octree/src/output.h"
 {
#if _MPI
    if (is_local(cell))
#endif
    {
#if _MPI
      if (fseek (p.fp, header + val(index,0,0,0)*cell_size, SEEK_SET) < 0) {
 perror ("output_gfs(): error while seeking");
 exit (1);
      }
#endif
      unsigned flags =
 level == 0 ? 0 :
#line 698 "/home/popinet/basilisk-octree/src/output.h"
      child.x == -1 && child.y == -1 && child.z == -1 ? 0 :
 child.x == -1 && child.y == -1 && child.z == 1 ? 1 :
 child.x == -1 && child.y == 1 && child.z == -1 ? 2 :
 child.x == -1 && child.y == 1 && child.z == 1 ? 3 :
 child.x == 1 && child.y == -1 && child.z == -1 ? 4 :
 child.x == 1 && child.y == -1 && child.z == 1 ? 5 :
 child.x == 1 && child.y == 1 && child.z == -1 ? 6 :
 7;

      if (is_leaf(cell))
 flags |= (1 << 4);
      fwrite (&flags, sizeof (unsigned), 1, p.fp);
      double a = -1;
      fwrite (&a, sizeof (double), 1, p.fp);
      if (list) for (scalar s = *list, *_i68 = list; ((scalar *)&s)->i >= 0; s = *++_i68)
 if (_attribute[s.i].name) {
   a = is_local(cell) && !isnan(val(s,0,0,0)) && val(s,0,0,0) != nodata ? val(s,0,0,0) : DBL_MAX;
   fwrite (&a, sizeof (double), 1, p.fp);
 }
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

#if _MPI
  delete (((scalar []){index,{-1}}));
  if (!pid() && fseek (p.fp, total_size, SEEK_SET) < 0) {
    perror ("output_gfs(): error while finishing");
    exit (1);
  }
  if (!pid())
#endif
    fputs ("}\n", p.fp);
  fflush (p.fp);

  if (!p.list)
    pfree (list,__func__,__FILE__,__LINE__);
  if (opened)
    fclose (p.fp);
 end_trace("output_gfs", "/home/popinet/basilisk-octree/src/output.h", 737); }
#line 760 "/home/popinet/basilisk-octree/src/output.h"
struct Dump {
  FILE * fp;
  scalar * list;
  double t;
  char * file;
};

struct DumpHeader {
  double t;
  long len;
  int depth;
};


void dump (struct Dump p)
{ trace ("dump", "/home/popinet/basilisk-octree/src/output.h", 775);
  FILE * fp = p.fp;
  scalar * lista = p.list ? p.list : all, * list = NULL;
  char * file = p.file;

  if (lista) for (scalar s = *lista, *_i69 = lista; ((scalar *)&s)->i >= 0; s = *++_i69)
    if (!_attribute[s.i].face && s.i != cm.i)
      list = list_add (list, s);

  if (file && (fp = fopen (file, "w")) == NULL) {
    perror (file);
    exit (1);
  }
  assert (fp);

  struct DumpHeader header = { p.t, list_len(list), depth() };

  if (pid() == 0 && fwrite (&header, sizeof(header), 1, fp) < 1) {
    perror ("dump(): error while writing header");
    exit (1);
  }

  scalar index = {-1};

#if _MPI
  index = new_scalar("index");
  z_indexing (index, false);
  int cell_size = sizeof(unsigned) + header.len*sizeof(double);
#endif

   { foreach_cell(){

#line 805 "/home/popinet/basilisk-octree/src/output.h"
 {
#if _MPI
    if (is_local(cell))
#endif
    {
#if _MPI
      if (fseek (fp, sizeof(header) + val(index,0,0,0)*cell_size, SEEK_SET) < 0) {
 perror ("dump(): error while seeking");
 exit (1);
      }
#endif
      unsigned flags = is_leaf(cell) ? leaf : 0;
      if (fwrite (&flags, sizeof(unsigned), 1, fp) < 1) {
 perror ("dump(): error while writing flags");
 exit (1);
      }
      if (list) for (scalar s = *list, *_i70 = list; ((scalar *)&s)->i >= 0; s = *++_i70)
 if (fwrite (&val(s,0,0,0), sizeof(double), 1, fp) < 1) {
   perror ("dump(): error while writing scalars");
   exit (1);
 }
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  delete (((scalar []){index,{-1}}));

  pfree (list,__func__,__FILE__,__LINE__);
  if (file)
    fclose (fp);
 end_trace("dump", "/home/popinet/basilisk-octree/src/output.h", 836); }


bool restore (struct Dump p)
{ trace ("restore", "/home/popinet/basilisk-octree/src/output.h", 840);
  FILE * fp = p.fp;
  scalar * lista = p.list ? p.list : all, * list = NULL;
  char * file = p.file;
  if (file && (fp = fopen (file, "r")) == NULL)
    { bool _ret =  false; end_trace("restore", "/home/popinet/basilisk-octree/src/output.h", 845);  return _ret; }
  assert (fp);

  if (lista) for (scalar s = *lista, *_i71 = lista; ((scalar *)&s)->i >= 0; s = *++_i71)
    if (!_attribute[s.i].face && s.i != cm.i)
      list = list_add (list, s);

  struct DumpHeader header;

  double t = 0.;
  if (fread (&header, sizeof(header), 1, fp) < 1) {
    fprintf (qstderr(), "restore(): error: expecting header\n");
    exit (1);
  }
  if (header.len != list_len (list)) {
    fprintf (qstderr(),
      "restore(): error: the list lengths don't match: %ld != %d\n",
      header.len, list_len (list));
    exit (1);
  }


  init_grid (1);

   { foreach_cell(){

#line 869 "/home/popinet/basilisk-octree/src/output.h"
 {
    cell.pid = pid();
    cell.flags |= active;
  } } end_foreach_cell(); }
  ((Quadtree *)grid)->dirty = true;




   { foreach_cell(){

#line 878 "/home/popinet/basilisk-octree/src/output.h"
 {
    unsigned flags;
    if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
      fprintf (qstderr(), "restore(): error: expecting 'flags'\n");
      exit (1);
    }
    if (list) for (scalar s = *list, *_i72 = list; ((scalar *)&s)->i >= 0; s = *++_i72)
      if (fread (&val(s,0,0,0), sizeof(double), 1, fp) != 1) {
 fprintf (qstderr(), "restore(): error: expecting '%s'\n", _attribute[s.i].name);
 exit (1);
      }
    if (!(flags & leaf) && is_leaf(cell))
      refine_cell (point, NULL, 0, NULL);
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }
  boundary (list);
  if (file)
    fclose (fp);


  double t1 = 0.;
  while (t1 < t && events (0, t1, false))
    t1 = tnext;

  pfree (list,__func__,__FILE__,__LINE__);
  { bool _ret =  true; end_trace("restore", "/home/popinet/basilisk-octree/src/output.h", 904);  return _ret; }
 end_trace("restore", "/home/popinet/basilisk-octree/src/output.h", 905); }
#line 277 "/home/popinet/basilisk-octree/src/utils.h"
#line 8 "/home/popinet/basilisk-octree/src/run.h"




double t = 0., dt = 1.;


void run (void)
{ trace ("run", "/home/popinet/basilisk-octree/src/run.h", 16);
  t = 0.; dt = 1.;
  init_grid (N);

  int i = 0;
  perf.nc = perf.tnc = 0;
  perf.gt = timer_start();
  while (events (i, t, true)) {





    update_perf();
    i++; t = tnext;
  }




  timer_print (perf.gt, i, perf.tnc);

  free_grid();
 end_trace("run", "/home/popinet/basilisk-octree/src/run.h", 39); }
#line 25 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
#line 1 "./timestep.h"
#line 1 "/home/popinet/basilisk-octree/src/timestep.h"

double timestep (const vector u, double dtmax)
{
  static double previous = 0.;
  dtmax /= CFL;
   { 
#undef _OMPSTART
#define _OMPSTART double _dtmax = dtmax; 
#undef _OMPEND
#define _OMPEND OMP(omp critical) if (_dtmax < dtmax) dtmax = _dtmax; mpi_all_reduce_double (dtmax, MPI_MIN); 
#line 6

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 6
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 6
{

#line 6 "/home/popinet/basilisk-octree/src/timestep.h"

    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta*val_cm(cm,0,0,0)/fabs(val(u.x,0,0,0));
      if (dt < _dtmax) _dtmax = dt;
    } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 6
{

#line 6 "/home/popinet/basilisk-octree/src/timestep.h"

    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta*val_cm(cm,0,0,0)/fabs(val(u.y,0,0,0));
      if (dt < _dtmax) _dtmax = dt;
    } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 6
{

#line 6 "/home/popinet/basilisk-octree/src/timestep.h"

    if (val(u.z,0,0,0) != 0.) {
      double dt = Delta*val_cm(cm,0,0,0)/fabs(val(u.z,0,0,0));
      if (dt < _dtmax) _dtmax = dt;
    } }  }}  end_foreach_face_generic()
#line 10
 end_foreach_face(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 6
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 6
{

#line 6 "/home/popinet/basilisk-octree/src/timestep.h"

    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta*val_cm(cm,0,0,0)/fabs(val(u.x,0,0,0));
      if (dt < _dtmax) _dtmax = dt;
    } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 6
{

#line 6 "/home/popinet/basilisk-octree/src/timestep.h"

    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta*val_cm(cm,0,0,0)/fabs(val(u.y,0,0,0));
      if (dt < _dtmax) _dtmax = dt;
    } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 6
{

#line 6 "/home/popinet/basilisk-octree/src/timestep.h"

    if (val(u.z,0,0,0) != 0.) {
      double dt = Delta*val_cm(cm,0,0,0)/fabs(val(u.z,0,0,0));
      if (dt < _dtmax) _dtmax = dt;
    } }  }}  end_foreach_face_generic()
#line 10
 end_foreach_face(); }
#undef _OMPSTART
#undef _OMPEND
#define _OMPSTART
#define _OMPEND
#line 11
 }
  dtmax *= CFL;
  if (dtmax > previous)
    dtmax = (previous + 0.1*dtmax)/1.1;
  previous = dtmax;
  return dtmax;
}
#line 26 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
#line 1 "./bcg.h"
#line 1 "/home/popinet/basilisk-octree/src/bcg.h"
#line 11 "/home/popinet/basilisk-octree/src/bcg.h"
void tracer_fluxes (scalar f,
      vector uf,
      vector flux,
      double dt,
       scalar src)
{





  vector g= new_vector("g");
  gradients (((scalar []){f,{-1}}), ((vector []){{g.x,g.y,g.z},{{-1},{-1},{-1}}}));




   { 
if (!is_constant(fm.x) && !is_constant(src)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_src
#define val_src(a,i,j,k) val(a,i,j,k)
#undef fine_src
#define fine_src(a,i,j,k) fine(a,i,j,k)
#undef coarse_src
#define coarse_src(a,i,j,k) coarse(a,i,j,k)
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 28
{

#line 28 "/home/popinet/basilisk-octree/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. +
      s*min(1., 1. - s*un)*val(g.x,i,0,0)*Delta/2.;





      double vn = val(uf.y,i,0,0)/val_fm_y(fm.y,i,0,0) + val(uf.y,i,1,0)/val_fm_y(fm.y,i,1,0);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.z,i,0,0)/val_fm_z(fm.z,i,0,0) + val(uf.z,i,0,1)/val_fm_z(fm.z,i,0,1);
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 28
{

#line 28 "/home/popinet/basilisk-octree/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. +
      s*min(1., 1. - s*un)*val(g.y,0,i,0)*Delta/2.;





      double vn = val(uf.z,0,i,0)/val_fm_z(fm.z,0,i,0) + val(uf.z,0,i,1)/val_fm_z(fm.z,0,i,1);
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.x,0,i,0)/val_fm_x(fm.x,0,i,0) + val(uf.x,1,i,0)/val_fm_x(fm.x,1,i,0);
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 28
{

#line 28 "/home/popinet/basilisk-octree/src/bcg.h"
 {







    double un = dt*val(uf.z,0,0,0)/(val_fm_z(fm.z,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (val_src(src,0,0,0) + val_src(src,0,0,-1))*dt/4. +
      s*min(1., 1. - s*un)*val(g.z,0,0,i)*Delta/2.;





      double vn = val(uf.x,0,0,i)/val_fm_x(fm.x,0,0,i) + val(uf.x,1,0,i)/val_fm_x(fm.x,1,0,i);
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.y,0,0,i)/val_fm_y(fm.y,0,0,i) + val(uf.y,0,1,i)/val_fm_y(fm.y,0,1,i);
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 56
 end_foreach_face(); }
if (is_constant(fm.x) && !is_constant(src)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_src
#define val_src(a,i,j,k) val(a,i,j,k)
#undef fine_src
#define fine_src(a,i,j,k) fine(a,i,j,k)
#undef coarse_src
#define coarse_src(a,i,j,k) coarse(a,i,j,k)
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 28
{

#line 28 "/home/popinet/basilisk-octree/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. +
      s*min(1., 1. - s*un)*val(g.x,i,0,0)*Delta/2.;





      double vn = val(uf.y,i,0,0)/val_fm_y(fm.y,i,0,0) + val(uf.y,i,1,0)/val_fm_y(fm.y,i,1,0);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.z,i,0,0)/val_fm_z(fm.z,i,0,0) + val(uf.z,i,0,1)/val_fm_z(fm.z,i,0,1);
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 28
{

#line 28 "/home/popinet/basilisk-octree/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. +
      s*min(1., 1. - s*un)*val(g.y,0,i,0)*Delta/2.;





      double vn = val(uf.z,0,i,0)/val_fm_z(fm.z,0,i,0) + val(uf.z,0,i,1)/val_fm_z(fm.z,0,i,1);
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.x,0,i,0)/val_fm_x(fm.x,0,i,0) + val(uf.x,1,i,0)/val_fm_x(fm.x,1,i,0);
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 28
{

#line 28 "/home/popinet/basilisk-octree/src/bcg.h"
 {







    double un = dt*val(uf.z,0,0,0)/(val_fm_z(fm.z,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (val_src(src,0,0,0) + val_src(src,0,0,-1))*dt/4. +
      s*min(1., 1. - s*un)*val(g.z,0,0,i)*Delta/2.;





      double vn = val(uf.x,0,0,i)/val_fm_x(fm.x,0,0,i) + val(uf.x,1,0,i)/val_fm_x(fm.x,1,0,i);
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.y,0,0,i)/val_fm_y(fm.y,0,0,i) + val(uf.y,0,1,i)/val_fm_y(fm.y,0,1,i);
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 56
 end_foreach_face(); }
if (!is_constant(fm.x) && is_constant(src)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const double _const_src = _constant[src.i -_NVARMAX];
NOT_UNUSED(_const_src);
#undef val_src
#define val_src(a,i,j,k) _const_src
#undef fine_src
#define fine_src(a,i,j,k) _const_src
#undef coarse_src
#define coarse_src(a,i,j,k) _const_src
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 28
{

#line 28 "/home/popinet/basilisk-octree/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. +
      s*min(1., 1. - s*un)*val(g.x,i,0,0)*Delta/2.;





      double vn = val(uf.y,i,0,0)/val_fm_y(fm.y,i,0,0) + val(uf.y,i,1,0)/val_fm_y(fm.y,i,1,0);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.z,i,0,0)/val_fm_z(fm.z,i,0,0) + val(uf.z,i,0,1)/val_fm_z(fm.z,i,0,1);
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 28
{

#line 28 "/home/popinet/basilisk-octree/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. +
      s*min(1., 1. - s*un)*val(g.y,0,i,0)*Delta/2.;





      double vn = val(uf.z,0,i,0)/val_fm_z(fm.z,0,i,0) + val(uf.z,0,i,1)/val_fm_z(fm.z,0,i,1);
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.x,0,i,0)/val_fm_x(fm.x,0,i,0) + val(uf.x,1,i,0)/val_fm_x(fm.x,1,i,0);
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 28
{

#line 28 "/home/popinet/basilisk-octree/src/bcg.h"
 {







    double un = dt*val(uf.z,0,0,0)/(val_fm_z(fm.z,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (val_src(src,0,0,0) + val_src(src,0,0,-1))*dt/4. +
      s*min(1., 1. - s*un)*val(g.z,0,0,i)*Delta/2.;





      double vn = val(uf.x,0,0,i)/val_fm_x(fm.x,0,0,i) + val(uf.x,1,0,i)/val_fm_x(fm.x,1,0,i);
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.y,0,0,i)/val_fm_y(fm.y,0,0,i) + val(uf.y,0,1,i)/val_fm_y(fm.y,0,1,i);
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 56
 end_foreach_face(); }
if (is_constant(fm.x) && is_constant(src)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const double _const_src = _constant[src.i -_NVARMAX];
NOT_UNUSED(_const_src);
#undef val_src
#define val_src(a,i,j,k) _const_src
#undef fine_src
#define fine_src(a,i,j,k) _const_src
#undef coarse_src
#define coarse_src(a,i,j,k) _const_src
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 28
{

#line 28 "/home/popinet/basilisk-octree/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. +
      s*min(1., 1. - s*un)*val(g.x,i,0,0)*Delta/2.;





      double vn = val(uf.y,i,0,0)/val_fm_y(fm.y,i,0,0) + val(uf.y,i,1,0)/val_fm_y(fm.y,i,1,0);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.z,i,0,0)/val_fm_z(fm.z,i,0,0) + val(uf.z,i,0,1)/val_fm_z(fm.z,i,0,1);
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 28
{

#line 28 "/home/popinet/basilisk-octree/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. +
      s*min(1., 1. - s*un)*val(g.y,0,i,0)*Delta/2.;





      double vn = val(uf.z,0,i,0)/val_fm_z(fm.z,0,i,0) + val(uf.z,0,i,1)/val_fm_z(fm.z,0,i,1);
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.x,0,i,0)/val_fm_x(fm.x,0,i,0) + val(uf.x,1,i,0)/val_fm_x(fm.x,1,i,0);
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 28
{

#line 28 "/home/popinet/basilisk-octree/src/bcg.h"
 {







    double un = dt*val(uf.z,0,0,0)/(val_fm_z(fm.z,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (val_src(src,0,0,0) + val_src(src,0,0,-1))*dt/4. +
      s*min(1., 1. - s*un)*val(g.z,0,0,i)*Delta/2.;





      double vn = val(uf.x,0,0,i)/val_fm_x(fm.x,0,0,i) + val(uf.x,1,0,i)/val_fm_x(fm.x,1,0,i);
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.y,0,0,i)/val_fm_y(fm.y,0,0,i) + val(uf.y,0,1,i)/val_fm_y(fm.y,0,1,i);
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 56
 end_foreach_face(); } }





  boundary_flux (((vector []){{flux.x,flux.y,flux.z},{{-1},{-1},{-1}}}));
 delete (((scalar []){g.x,g.y,g.z,{-1}})); }






struct Advection {
  scalar * tracers;
  vector u;
  double dt;
  scalar * src;
};

void advection (struct Advection p)
{




  scalar * lsrc = p.src;
  if (!lsrc) {
    scalar zero= new_const_scalar("zero", 8,  0.);
    if (p.tracers) for (scalar s = *p.tracers, *_i73 = p.tracers; ((scalar *)&s)->i >= 0; s = *++_i73)
      lsrc = list_append (lsrc, zero);
  }

  assert (list_len(p.tracers) == list_len(lsrc));
  scalar f, src;
  scalar * _i2 = p.tracers; scalar * _i3 = lsrc; if (p.tracers) for (f = *p.tracers, src = *lsrc; ((scalar *)&f)->i >= 0; f = *++_i2, src = *++_i3) {
    vector flux= new_face_vector("flux");
    tracer_fluxes (f, p.u, flux, p.dt, src);
     { 
if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 95
foreach(){

#line 95 "/home/popinet/basilisk-octree/src/bcg.h"

      {
#line 96

        val(f,0,0,0) += p.dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*val_cm(cm,0,0,0));
#line 96

        val(f,0,0,0) += p.dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*val_cm(cm,0,0,0));
#line 96

        val(f,0,0,0) += p.dt*(val(flux.z,0,0,0) - val(flux.z,0,0,1))/(Delta*val_cm(cm,0,0,0));}; } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 95
foreach(){

#line 95 "/home/popinet/basilisk-octree/src/bcg.h"

      {
#line 96

        val(f,0,0,0) += p.dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*val_cm(cm,0,0,0));
#line 96

        val(f,0,0,0) += p.dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*val_cm(cm,0,0,0));
#line 96

        val(f,0,0,0) += p.dt*(val(flux.z,0,0,0) - val(flux.z,0,0,1))/(Delta*val_cm(cm,0,0,0));}; } end_foreach(); } }
    boundary (((scalar []){f,{-1}}));
   delete (((scalar []){flux.x,flux.y,flux.z,{-1}})); }

  if (!p.src)
    pfree (lsrc,__func__,__FILE__,__LINE__);
}
#line 27 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
#line 1 "./viscosity.h"
#line 1 "/home/popinet/basilisk-octree/src/viscosity.h"
#line 1 "./poisson.h"
#line 1 "/home/popinet/basilisk-octree/src/poisson.h"
#line 32 "/home/popinet/basilisk-octree/src/poisson.h"
void mg_cycle (scalar * a, scalar * res, scalar * da,
        void (* relax) (scalar * da, scalar * res,
          int depth, void * data),
        void * data,
        int nrelax, int minlevel, int maxlevel)
{




  restriction (res);





  for (int l = minlevel; l <= maxlevel; l++) {




    if (l == minlevel)
       { foreach_level_or_leaf (l){

#line 54 "/home/popinet/basilisk-octree/src/poisson.h"

 if (da) for (scalar s = *da, *_i74 = da; ((scalar *)&s)->i >= 0; s = *++_i74)
   val(s,0,0,0) = 0.; } end_foreach_level_or_leaf(); }





    else
       { foreach_level (l){

#line 63 "/home/popinet/basilisk-octree/src/poisson.h"

 if (da) for (scalar s = *da, *_i75 = da; ((scalar *)&s)->i >= 0; s = *++_i75)
   val(s,0,0,0) = bilinear (point, s); } end_foreach_level(); }





    boundary_level (da, l);
    for (int i = 0; i < nrelax; i++) {
      relax (da, res, l, data);
      boundary_level (da, l);
    }
  }




   { foreach(){

#line 81 "/home/popinet/basilisk-octree/src/poisson.h"
 {
    scalar s, ds;
    scalar * _i4 = a; scalar * _i5 = da; if (a) for (s = *a, ds = *da; ((scalar *)&s)->i >= 0; s = *++_i4, ds = *++_i5)
      val(s,0,0,0) += val(ds,0,0,0);
  } } end_foreach(); }
  boundary (a);
}
#line 99 "/home/popinet/basilisk-octree/src/poisson.h"
int NITERMAX = 100;
double TOLERANCE = 1e-3;




typedef struct {
  int i;
  double resb, resa;
  double sum;
} mgstats;







mgstats mg_solve (scalar * a, scalar * b,
    double (* residual) (scalar * a, scalar * b, scalar * res,
           void * data),
    void (* relax) (scalar * da, scalar * res, int depth,
      void * data),
    void * data)
{





  scalar * da = list_clone (a), * res = NULL;
  if (a) for (scalar s = *a, *_i76 = a; ((scalar *)&s)->i >= 0; s = *++_i76) {
    scalar r = new_scalar("r");
    res = list_append (res, r);
  }






  for (int b = 0; b < nboundary; b++)
    if (da) for (scalar s = *da, *_i77 = da; ((scalar *)&s)->i >= 0; s = *++_i77)
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b];




  mgstats s = {0, 0., 0.};
  double sum = 0.;
   { 
#undef _OMPSTART
#define _OMPSTART 
#undef _OMPEND
#define _OMPEND mpi_all_reduce_double (sum, MPI_SUM); 
#line 149
foreach (reduction(+:sum)){

#line 149 "/home/popinet/basilisk-octree/src/poisson.h"

    if (b) for (scalar s = *b, *_i78 = b; ((scalar *)&s)->i >= 0; s = *++_i78)
      sum += val(s,0,0,0); } end_foreach();
#undef _OMPSTART
#undef _OMPEND
#define _OMPSTART
#define _OMPEND
#line 152
 }
  s.sum = sum;




  s.resb = s.resa = residual (a, b, res, data);




  int maxlevel = depth();
  mpi_all_reduce (maxlevel, MPI_INT, MPI_MAX);






  for (s.i = 0; s.i < NITERMAX && (s.i < 1 || s.resa > TOLERANCE); s.i++) {
    mg_cycle (a, res, da, relax, data, 4, 0, maxlevel);
    s.resa = residual (a, b, res, data);
  }




  if (s.i == NITERMAX)
    fprintf (ferr,
      "WARNING: convergence not reached after %d iterations\n"
      "  res: %g sum: %g\n",
      NITERMAX, s.resa, s.sum), fflush (ferr);




  delete (res); pfree (res,__func__,__FILE__,__LINE__);
  delete (da); pfree (da,__func__,__FILE__,__LINE__);

  return s;
}
#line 210 "/home/popinet/basilisk-octree/src/poisson.h"
struct Poisson {
  scalar a, b;
   vector alpha;
   scalar lambda;
  double tolerance;
};





static void relax (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson * p = data;
   vector alpha = p->alpha;
   scalar lambda = p->lambda;
#line 243 "/home/popinet/basilisk-octree/src/poisson.h"
  scalar c = a;






   { 
if (!is_constant(lambda) && !is_constant(alpha.x)) {
#undef val_lambda
#define val_lambda(a,i,j,k) val(a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 250
foreach_level_or_leaf (l){

#line 250 "/home/popinet/basilisk-octree/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 252
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 252
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }
#line 252
 {
      n += val_alpha_z(alpha.z,0,0,1)*val(a,0,0,1) + val_alpha_z(alpha.z,0,0,0)*val(a,0,0,-1);
      d += val_alpha_z(alpha.z,0,0,1) + val_alpha_z(alpha.z,0,0,0);
    }}
    val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
if (is_constant(lambda) && !is_constant(alpha.x)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 250
foreach_level_or_leaf (l){

#line 250 "/home/popinet/basilisk-octree/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 252
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 252
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }
#line 252
 {
      n += val_alpha_z(alpha.z,0,0,1)*val(a,0,0,1) + val_alpha_z(alpha.z,0,0,0)*val(a,0,0,-1);
      d += val_alpha_z(alpha.z,0,0,1) + val_alpha_z(alpha.z,0,0,0);
    }}
    val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
if (!is_constant(lambda) && is_constant(alpha.x)) {
#undef val_lambda
#define val_lambda(a,i,j,k) val(a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 250
foreach_level_or_leaf (l){

#line 250 "/home/popinet/basilisk-octree/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 252
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 252
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }
#line 252
 {
      n += val_alpha_z(alpha.z,0,0,1)*val(a,0,0,1) + val_alpha_z(alpha.z,0,0,0)*val(a,0,0,-1);
      d += val_alpha_z(alpha.z,0,0,1) + val_alpha_z(alpha.z,0,0,0);
    }}
    val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
if (is_constant(lambda) && is_constant(alpha.x)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 250
foreach_level_or_leaf (l){

#line 250 "/home/popinet/basilisk-octree/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 252
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 252
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }
#line 252
 {
      n += val_alpha_z(alpha.z,0,0,1)*val(a,0,0,1) + val_alpha_z(alpha.z,0,0,0)*val(a,0,0,-1);
      d += val_alpha_z(alpha.z,0,0,1) + val_alpha_z(alpha.z,0,0,0);
    }}
    val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); } }
#line 275 "/home/popinet/basilisk-octree/src/poisson.h"
}






static double residual (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson * p = data;
   vector alpha = p->alpha;
   scalar lambda = p->lambda;
  double maxres = 0.;


  vector g= new_face_vector("g");
   { 
if (!is_constant(alpha.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 292
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 292
{

#line 292 "/home/popinet/basilisk-octree/src/poisson.h"

    val(g.x,0,0,0) = val_alpha_x(alpha.x,0,0,0)*(val(a,0,0,0) - val(a,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 292
{

#line 292 "/home/popinet/basilisk-octree/src/poisson.h"

    val(g.y,0,0,0) = val_alpha_y(alpha.y,0,0,0)*(val(a,0,0,0) - val(a,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 292
{

#line 292 "/home/popinet/basilisk-octree/src/poisson.h"

    val(g.z,0,0,0) = val_alpha_z(alpha.z,0,0,0)*(val(a,0,0,0) - val(a,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 293
 end_foreach_face(); }
if (is_constant(alpha.x)) {
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 292
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 292
{

#line 292 "/home/popinet/basilisk-octree/src/poisson.h"

    val(g.x,0,0,0) = val_alpha_x(alpha.x,0,0,0)*(val(a,0,0,0) - val(a,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 292
{

#line 292 "/home/popinet/basilisk-octree/src/poisson.h"

    val(g.y,0,0,0) = val_alpha_y(alpha.y,0,0,0)*(val(a,0,0,0) - val(a,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 292
{

#line 292 "/home/popinet/basilisk-octree/src/poisson.h"

    val(g.z,0,0,0) = val_alpha_z(alpha.z,0,0,0)*(val(a,0,0,0) - val(a,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 293
 end_foreach_face(); } }
  boundary_flux (((vector []){{g.x,g.y,g.z},{{-1},{-1},{-1}}}));
   { 
#undef _OMPSTART
#define _OMPSTART double _maxres = maxres; 
#undef _OMPEND
#define _OMPEND OMP(omp critical) if (_maxres > maxres) maxres = _maxres; mpi_all_reduce_double (maxres, MPI_MAX); 
#line 295

if (!is_constant(lambda)) {
#undef val_lambda
#define val_lambda(a,i,j,k) val(a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) coarse(a,i,j,k)
#line 295
foreach (){

#line 295 "/home/popinet/basilisk-octree/src/poisson.h"
 {
    val(res,0,0,0) = val(b,0,0,0) - val_lambda(lambda,0,0,0)*val(a,0,0,0);
    {
#line 297

      val(res,0,0,0) += (val(g.x,0,0,0) - val(g.x,1,0,0))/Delta;
#line 297

      val(res,0,0,0) += (val(g.y,0,0,0) - val(g.y,0,1,0))/Delta;
#line 297

      val(res,0,0,0) += (val(g.z,0,0,0) - val(g.z,0,0,1))/Delta;}
    if (fabs (val(res,0,0,0)) > _maxres)
      _maxres = fabs (val(res,0,0,0));
  } } end_foreach(); }
if (is_constant(lambda)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
#line 295
foreach (){

#line 295 "/home/popinet/basilisk-octree/src/poisson.h"
 {
    val(res,0,0,0) = val(b,0,0,0) - val_lambda(lambda,0,0,0)*val(a,0,0,0);
    {
#line 297

      val(res,0,0,0) += (val(g.x,0,0,0) - val(g.x,1,0,0))/Delta;
#line 297

      val(res,0,0,0) += (val(g.y,0,0,0) - val(g.y,0,1,0))/Delta;
#line 297

      val(res,0,0,0) += (val(g.z,0,0,0) - val(g.z,0,0,1))/Delta;}
    if (fabs (val(res,0,0,0)) > _maxres)
      _maxres = fabs (val(res,0,0,0));
  } } end_foreach(); }
#undef _OMPSTART
#undef _OMPEND
#define _OMPSTART
#define _OMPEND
#line 302
 }
#line 313 "/home/popinet/basilisk-octree/src/poisson.h"
  boundary (resl);
  { double _ret =  maxres; delete (((scalar []){g.x,g.y,g.z,{-1}}));  return _ret; }
 delete (((scalar []){g.x,g.y,g.z,{-1}})); }
#line 326 "/home/popinet/basilisk-octree/src/poisson.h"
mgstats poisson (struct Poisson p)
{






  if (!p.alpha.x.i) {
    vector alpha= new_const_vector("alpha", 9, (double []) {1.,1.,1.});
    p.alpha = alpha;
  }
  if (!p.lambda.i) {
    scalar lambda= new_const_scalar("lambda", 12,  0.);
    p.lambda = lambda;
  }




  vector alpha = p.alpha;
  scalar lambda = p.lambda;
  restriction (((scalar []){alpha.x,alpha.y,alpha.z,lambda,{-1}}));





  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;

  scalar a = p.a, b = p.b;
  mgstats s = mg_solve (((scalar []){a,{-1}}), ((scalar []){b,{-1}}), residual, relax, &p);




  if (p.tolerance)
    TOLERANCE = defaultol;

  return s;
}
#line 387 "/home/popinet/basilisk-octree/src/poisson.h"

mgstats project (vector u, scalar p,  vector alpha, double dt)
{ trace ("project", "/home/popinet/basilisk-octree/src/poisson.h", 389);






  scalar div= new_scalar("div");
   { foreach(){

#line 397 "/home/popinet/basilisk-octree/src/poisson.h"
 {
    val(div,0,0,0) = 0.;
    {
#line 399

      val(div,0,0,0) += val(u.x,1,0,0) - val(u.x,0,0,0);
#line 399

      val(div,0,0,0) += val(u.y,0,1,0) - val(u.y,0,0,0);
#line 399

      val(div,0,0,0) += val(u.z,0,0,1) - val(u.z,0,0,0);}
    val(div,0,0,0) /= dt*Delta;
  } } end_foreach(); }
#line 413 "/home/popinet/basilisk-octree/src/poisson.h"
  mgstats mgp = poisson ((struct Poisson){p, div, alpha, .tolerance = TOLERANCE/sq(dt)});




   { 
if (!is_constant(alpha.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 418
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 418
{

#line 418 "/home/popinet/basilisk-octree/src/poisson.h"

    val(u.x,0,0,0) -= dt*val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 418
{

#line 418 "/home/popinet/basilisk-octree/src/poisson.h"

    val(u.y,0,0,0) -= dt*val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 418
{

#line 418 "/home/popinet/basilisk-octree/src/poisson.h"

    val(u.z,0,0,0) -= dt*val_alpha_z(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 419
 end_foreach_face(); }
if (is_constant(alpha.x)) {
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 418
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 418
{

#line 418 "/home/popinet/basilisk-octree/src/poisson.h"

    val(u.x,0,0,0) -= dt*val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 418
{

#line 418 "/home/popinet/basilisk-octree/src/poisson.h"

    val(u.y,0,0,0) -= dt*val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 418
{

#line 418 "/home/popinet/basilisk-octree/src/poisson.h"

    val(u.z,0,0,0) -= dt*val_alpha_z(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 419
 end_foreach_face(); } }
  boundary ((scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}));

  { mgstats _ret =  mgp; delete (((scalar []){div,{-1}}));  end_trace("project", "/home/popinet/basilisk-octree/src/poisson.h", 422);  return _ret; }
 delete (((scalar []){div,{-1}}));  end_trace("project", "/home/popinet/basilisk-octree/src/poisson.h", 423); }
#line 2 "/home/popinet/basilisk-octree/src/viscosity.h"

struct Viscosity {
  vector u;
  vector mu;
  scalar rho;
  double dt;
};
#line 23 "/home/popinet/basilisk-octree/src/viscosity.h"
static void relax_viscosity (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = data;
   vector mu = p->mu;
   scalar rho = p->rho;
  double dt = p->dt;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0])));




  vector w = u;


   { 
if (!is_constant(rho) && !is_constant(mu.x)) {
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_z
#define val_mu_z(a,i,j,k) val(a,i,j,k)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,i,j,k)
#line 37
foreach_level_or_leaf (l){

#line 37 "/home/popinet/basilisk-octree/src/viscosity.h"
 {
    {
#line 38

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


      + val_mu_z(mu.z,0,0,1)*(val(u.x,0,0,1) +
       (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
       (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.x,0,0,-1) +
         (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
         (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)


          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)

        ));
#line 38

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_z(mu.z,0,0,1)*(val(u.y,0,0,1) +
     (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
     (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.y,0,0,-1) +
         (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
         (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
       (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
       (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)


          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)

        ));
#line 38

      val(w.z,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1)*val(u.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)*val(u.z,0,0,-1)

      + val_mu_x(mu.x,1,0,0)*(val(u.z,1,0,0) +
     (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
     (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.z,-1,0,0) +
         (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
         (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


      + val_mu_y(mu.y,0,1,0)*(val(u.z,0,1,0) +
       (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
       (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.z,0,-1,0) +
         (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
         (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

      ) + val(r.z,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).z + dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)


          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)

        ));}
  } } end_foreach_level_or_leaf(); }
if (is_constant(rho) && !is_constant(mu.x)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_z
#define val_mu_z(a,i,j,k) val(a,i,j,k)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,i,j,k)
#line 37
foreach_level_or_leaf (l){

#line 37 "/home/popinet/basilisk-octree/src/viscosity.h"
 {
    {
#line 38

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


      + val_mu_z(mu.z,0,0,1)*(val(u.x,0,0,1) +
       (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
       (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.x,0,0,-1) +
         (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
         (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)


          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)

        ));
#line 38

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_z(mu.z,0,0,1)*(val(u.y,0,0,1) +
     (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
     (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.y,0,0,-1) +
         (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
         (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
       (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
       (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)


          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)

        ));
#line 38

      val(w.z,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1)*val(u.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)*val(u.z,0,0,-1)

      + val_mu_x(mu.x,1,0,0)*(val(u.z,1,0,0) +
     (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
     (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.z,-1,0,0) +
         (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
         (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


      + val_mu_y(mu.y,0,1,0)*(val(u.z,0,1,0) +
       (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
       (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.z,0,-1,0) +
         (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
         (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

      ) + val(r.z,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).z + dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)


          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)

        ));}
  } } end_foreach_level_or_leaf(); }
if (!is_constant(rho) && is_constant(mu.x)) {
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX], _constant[mu.z.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,i,j,k) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#line 37
foreach_level_or_leaf (l){

#line 37 "/home/popinet/basilisk-octree/src/viscosity.h"
 {
    {
#line 38

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


      + val_mu_z(mu.z,0,0,1)*(val(u.x,0,0,1) +
       (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
       (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.x,0,0,-1) +
         (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
         (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)


          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)

        ));
#line 38

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_z(mu.z,0,0,1)*(val(u.y,0,0,1) +
     (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
     (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.y,0,0,-1) +
         (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
         (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
       (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
       (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)


          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)

        ));
#line 38

      val(w.z,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1)*val(u.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)*val(u.z,0,0,-1)

      + val_mu_x(mu.x,1,0,0)*(val(u.z,1,0,0) +
     (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
     (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.z,-1,0,0) +
         (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
         (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


      + val_mu_y(mu.y,0,1,0)*(val(u.z,0,1,0) +
       (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
       (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.z,0,-1,0) +
         (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
         (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

      ) + val(r.z,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).z + dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)


          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)

        ));}
  } } end_foreach_level_or_leaf(); }
if (is_constant(rho) && is_constant(mu.x)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
const struct { double x, y, z; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX], _constant[mu.z.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,i,j,k) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#line 37
foreach_level_or_leaf (l){

#line 37 "/home/popinet/basilisk-octree/src/viscosity.h"
 {
    {
#line 38

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


      + val_mu_z(mu.z,0,0,1)*(val(u.x,0,0,1) +
       (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
       (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.x,0,0,-1) +
         (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
         (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)


          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)

        ));
#line 38

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_z(mu.z,0,0,1)*(val(u.y,0,0,1) +
     (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
     (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.y,0,0,-1) +
         (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
         (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
       (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
       (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)


          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)

        ));
#line 38

      val(w.z,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1)*val(u.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)*val(u.z,0,0,-1)

      + val_mu_x(mu.x,1,0,0)*(val(u.z,1,0,0) +
     (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
     (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.z,-1,0,0) +
         (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
         (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


      + val_mu_y(mu.y,0,1,0)*(val(u.z,0,1,0) +
       (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
       (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.z,0,-1,0) +
         (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
         (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

      ) + val(r.z,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).z + dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)


          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)

        ));}
  } } end_foreach_level_or_leaf(); } }
#line 83 "/home/popinet/basilisk-octree/src/viscosity.h"
}

static double residual_viscosity (scalar * a, scalar * b, scalar * resl,
      void * data)
{
  struct Viscosity * p = data;
   vector mu = p->mu;
   scalar rho = p->rho;
  double dt = p->dt;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0]))), res = (*((vector *)&(resl[0])));
  double maxres = 0.;


  {
#line 96
 {
    vector taux= new_face_vector("taux");
     { 
if (!is_constant(mu.x)) {
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_z
#define val_mu_z(a,i,j,k) val(a,i,j,k)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,i,j,k)
#line 98
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 98
{

#line 98 "/home/popinet/basilisk-octree/src/viscosity.h"

      val(taux.x,0,0,0) = 2.*val_mu_x(mu.x,0,0,0)*(val(u.x,0,0,0) - val(u.x,-1,0,0))/Delta; }  }}  end_foreach_face_generic()
#line 99
 end_foreach_face(); }
if (is_constant(mu.x)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX], _constant[mu.z.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,i,j,k) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#line 98
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 98
{

#line 98 "/home/popinet/basilisk-octree/src/viscosity.h"

      val(taux.x,0,0,0) = 2.*val_mu_x(mu.x,0,0,0)*(val(u.x,0,0,0) - val(u.x,-1,0,0))/Delta; }  }}  end_foreach_face_generic()
#line 99
 end_foreach_face(); } }

       { 
if (!is_constant(mu.x)) {
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_z
#define val_mu_z(a,i,j,k) val(a,i,j,k)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,i,j,k)
#line 101
foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 101
{

#line 101 "/home/popinet/basilisk-octree/src/viscosity.h"

 val(taux.y,0,0,0) = val_mu_y(mu.y,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
      (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
      (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 104
 end_foreach_face(); }
if (is_constant(mu.x)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX], _constant[mu.z.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,i,j,k) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#line 101
foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 101
{

#line 101 "/home/popinet/basilisk-octree/src/viscosity.h"

 val(taux.y,0,0,0) = val_mu_y(mu.y,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
      (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
      (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 104
 end_foreach_face(); } }


       { 
if (!is_constant(mu.x)) {
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_z
#define val_mu_z(a,i,j,k) val(a,i,j,k)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,i,j,k)
#line 107
foreach_face_generic() { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 107
{

#line 107 "/home/popinet/basilisk-octree/src/viscosity.h"

 val(taux.z,0,0,0) = val_mu_z(mu.z,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,0,-1) +
      (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
      (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 110
 end_foreach_face(); }
if (is_constant(mu.x)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX], _constant[mu.z.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,i,j,k) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#line 107
foreach_face_generic() { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 107
{

#line 107 "/home/popinet/basilisk-octree/src/viscosity.h"

 val(taux.z,0,0,0) = val_mu_z(mu.z,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,0,-1) +
      (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
      (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 110
 end_foreach_face(); } }

    boundary_flux (((vector []){{taux.x,taux.y,taux.z},{{-1},{-1},{-1}}}));
     { 
#undef _OMPSTART
#define _OMPSTART double _maxres = maxres; 
#undef _OMPEND
#define _OMPEND OMP(omp critical) if (_maxres > maxres) maxres = _maxres; mpi_all_reduce_double (maxres, MPI_MAX); 
#line 113

if (!is_constant(rho)) {
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
#line 113
foreach (){

#line 113 "/home/popinet/basilisk-octree/src/viscosity.h"
 {

      double d = 0.;
      d += val(taux.x,1,0,0) - val(taux.x,0,0,0);

        d += val(taux.y,0,1,0) - val(taux.y,0,0,0);


        d += val(taux.z,0,0,1) - val(taux.z,0,0,0);

      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.,1.}).x*val(u.x,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.x,0,0,0)) > _maxres)
 _maxres = fabs (val(res.x,0,0,0));
    } } end_foreach(); }
if (is_constant(rho)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#line 113
foreach (){

#line 113 "/home/popinet/basilisk-octree/src/viscosity.h"
 {

      double d = 0.;
      d += val(taux.x,1,0,0) - val(taux.x,0,0,0);

        d += val(taux.y,0,1,0) - val(taux.y,0,0,0);


        d += val(taux.z,0,0,1) - val(taux.z,0,0,0);

      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.,1.}).x*val(u.x,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.x,0,0,0)) > _maxres)
 _maxres = fabs (val(res.x,0,0,0));
    } } end_foreach(); }
#undef _OMPSTART
#undef _OMPEND
#define _OMPSTART
#define _OMPEND
#line 127
 }
   delete (((scalar []){taux.x,taux.y,taux.z,{-1}})); }
#line 96
 {
    vector taux= new_face_vector("taux");
     { 
if (!is_constant(mu.y)) {
#undef val_mu_y
#define val_mu_y(a,k,i,j) val(a,k,i,j)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,k,i,j)
#undef val_mu_z
#define val_mu_z(a,k,i,j) val(a,k,i,j)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,k,i,j)
#undef val_mu_x
#define val_mu_x(a,k,i,j) val(a,k,i,j)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,k,i,j)
#line 98
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_y()) {
#line 98
{

#line 98 "/home/popinet/basilisk-octree/src/viscosity.h"

      val(taux.y,0,0,0) = 2.*val_mu_y(mu.y,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 99
 end_foreach_face(); }
if (is_constant(mu.y)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.y.i -_NVARMAX], _constant[mu.z.i - _NVARMAX], _constant[mu.x.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_y
#define val_mu_y(a,k,i,j) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,k,i,j) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#undef val_mu_x
#define val_mu_x(a,k,i,j) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#line 98
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_y()) {
#line 98
{

#line 98 "/home/popinet/basilisk-octree/src/viscosity.h"

      val(taux.y,0,0,0) = 2.*val_mu_y(mu.y,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 99
 end_foreach_face(); } }

       { 
if (!is_constant(mu.y)) {
#undef val_mu_y
#define val_mu_y(a,k,i,j) val(a,k,i,j)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,k,i,j)
#undef val_mu_z
#define val_mu_z(a,k,i,j) val(a,k,i,j)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,k,i,j)
#undef val_mu_x
#define val_mu_x(a,k,i,j) val(a,k,i,j)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,k,i,j)
#line 101
foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_z()) {
#line 101
{

#line 101 "/home/popinet/basilisk-octree/src/viscosity.h"

 val(taux.z,0,0,0) = val_mu_z(mu.z,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,0,-1) +
      (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
      (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 104
 end_foreach_face(); }
if (is_constant(mu.y)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.y.i -_NVARMAX], _constant[mu.z.i - _NVARMAX], _constant[mu.x.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_y
#define val_mu_y(a,k,i,j) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,k,i,j) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#undef val_mu_x
#define val_mu_x(a,k,i,j) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#line 101
foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_z()) {
#line 101
{

#line 101 "/home/popinet/basilisk-octree/src/viscosity.h"

 val(taux.z,0,0,0) = val_mu_z(mu.z,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,0,-1) +
      (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
      (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 104
 end_foreach_face(); } }


       { 
if (!is_constant(mu.y)) {
#undef val_mu_y
#define val_mu_y(a,k,i,j) val(a,k,i,j)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,k,i,j)
#undef val_mu_z
#define val_mu_z(a,k,i,j) val(a,k,i,j)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,k,i,j)
#undef val_mu_x
#define val_mu_x(a,k,i,j) val(a,k,i,j)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,k,i,j)
#line 107
foreach_face_generic() { int kg = -1; VARIABLES;  if (is_face_x()) {
#line 107
{

#line 107 "/home/popinet/basilisk-octree/src/viscosity.h"

 val(taux.x,0,0,0) = val_mu_x(mu.x,0,0,0)*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
      (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
      (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 110
 end_foreach_face(); }
if (is_constant(mu.y)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.y.i -_NVARMAX], _constant[mu.z.i - _NVARMAX], _constant[mu.x.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_y
#define val_mu_y(a,k,i,j) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,k,i,j) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#undef val_mu_x
#define val_mu_x(a,k,i,j) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#line 107
foreach_face_generic() { int kg = -1; VARIABLES;  if (is_face_x()) {
#line 107
{

#line 107 "/home/popinet/basilisk-octree/src/viscosity.h"

 val(taux.x,0,0,0) = val_mu_x(mu.x,0,0,0)*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
      (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
      (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 110
 end_foreach_face(); } }

    boundary_flux (((vector []){{taux.x,taux.y,taux.z},{{-1},{-1},{-1}}}));
     { 
#undef _OMPSTART
#define _OMPSTART double _maxres = maxres; 
#undef _OMPEND
#define _OMPEND OMP(omp critical) if (_maxres > maxres) maxres = _maxres; mpi_all_reduce_double (maxres, MPI_MAX); 
#line 113

if (!is_constant(rho)) {
#undef val_rho
#define val_rho(a,k,i,j) val(a,k,i,j)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,k,i,j)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,k,i,j)
#line 113
foreach (){

#line 113 "/home/popinet/basilisk-octree/src/viscosity.h"
 {

      double d = 0.;
      d += val(taux.y,0,1,0) - val(taux.y,0,0,0);

        d += val(taux.z,0,0,1) - val(taux.z,0,0,0);


        d += val(taux.x,1,0,0) - val(taux.x,0,0,0);

      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.,1.}).y*val(u.y,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.y,0,0,0)) > _maxres)
 _maxres = fabs (val(res.y,0,0,0));
    } } end_foreach(); }
if (is_constant(rho)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,k,i,j) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#line 113
foreach (){

#line 113 "/home/popinet/basilisk-octree/src/viscosity.h"
 {

      double d = 0.;
      d += val(taux.y,0,1,0) - val(taux.y,0,0,0);

        d += val(taux.z,0,0,1) - val(taux.z,0,0,0);


        d += val(taux.x,1,0,0) - val(taux.x,0,0,0);

      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.,1.}).y*val(u.y,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.y,0,0,0)) > _maxres)
 _maxres = fabs (val(res.y,0,0,0));
    } } end_foreach(); }
#undef _OMPSTART
#undef _OMPEND
#define _OMPSTART
#define _OMPEND
#line 127
 }
   delete (((scalar []){taux.x,taux.y,taux.z,{-1}})); }
#line 96
 {
    vector taux= new_face_vector("taux");
     { 
if (!is_constant(mu.z)) {
#undef val_mu_z
#define val_mu_z(a,j,k,i) val(a,j,k,i)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,j,k,i)
#undef val_mu_x
#define val_mu_x(a,j,k,i) val(a,j,k,i)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,j,k,i)
#undef val_mu_y
#define val_mu_y(a,j,k,i) val(a,j,k,i)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,j,k,i)
#line 98
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_z()) {
#line 98
{

#line 98 "/home/popinet/basilisk-octree/src/viscosity.h"

      val(taux.z,0,0,0) = 2.*val_mu_z(mu.z,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 99
 end_foreach_face(); }
if (is_constant(mu.z)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.z.i -_NVARMAX], _constant[mu.x.i - _NVARMAX], _constant[mu.y.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_z
#define val_mu_z(a,j,k,i) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#undef val_mu_x
#define val_mu_x(a,j,k,i) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,j,k,i) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#line 98
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_z()) {
#line 98
{

#line 98 "/home/popinet/basilisk-octree/src/viscosity.h"

      val(taux.z,0,0,0) = 2.*val_mu_z(mu.z,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 99
 end_foreach_face(); } }

       { 
if (!is_constant(mu.z)) {
#undef val_mu_z
#define val_mu_z(a,j,k,i) val(a,j,k,i)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,j,k,i)
#undef val_mu_x
#define val_mu_x(a,j,k,i) val(a,j,k,i)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,j,k,i)
#undef val_mu_y
#define val_mu_y(a,j,k,i) val(a,j,k,i)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,j,k,i)
#line 101
foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_x()) {
#line 101
{

#line 101 "/home/popinet/basilisk-octree/src/viscosity.h"

 val(taux.x,0,0,0) = val_mu_x(mu.x,0,0,0)*(val(u.z,0,0,0) - val(u.z,-1,0,0) +
      (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
      (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 104
 end_foreach_face(); }
if (is_constant(mu.z)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.z.i -_NVARMAX], _constant[mu.x.i - _NVARMAX], _constant[mu.y.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_z
#define val_mu_z(a,j,k,i) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#undef val_mu_x
#define val_mu_x(a,j,k,i) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,j,k,i) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#line 101
foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_x()) {
#line 101
{

#line 101 "/home/popinet/basilisk-octree/src/viscosity.h"

 val(taux.x,0,0,0) = val_mu_x(mu.x,0,0,0)*(val(u.z,0,0,0) - val(u.z,-1,0,0) +
      (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
      (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 104
 end_foreach_face(); } }


       { 
if (!is_constant(mu.z)) {
#undef val_mu_z
#define val_mu_z(a,j,k,i) val(a,j,k,i)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,j,k,i)
#undef val_mu_x
#define val_mu_x(a,j,k,i) val(a,j,k,i)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,j,k,i)
#undef val_mu_y
#define val_mu_y(a,j,k,i) val(a,j,k,i)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,j,k,i)
#line 107
foreach_face_generic() { int kg = -1; VARIABLES;  if (is_face_y()) {
#line 107
{

#line 107 "/home/popinet/basilisk-octree/src/viscosity.h"

 val(taux.y,0,0,0) = val_mu_y(mu.y,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,-1,0) +
      (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
      (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 110
 end_foreach_face(); }
if (is_constant(mu.z)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.z.i -_NVARMAX], _constant[mu.x.i - _NVARMAX], _constant[mu.y.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_z
#define val_mu_z(a,j,k,i) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#undef val_mu_x
#define val_mu_x(a,j,k,i) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,j,k,i) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#line 107
foreach_face_generic() { int kg = -1; VARIABLES;  if (is_face_y()) {
#line 107
{

#line 107 "/home/popinet/basilisk-octree/src/viscosity.h"

 val(taux.y,0,0,0) = val_mu_y(mu.y,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,-1,0) +
      (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
      (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 110
 end_foreach_face(); } }

    boundary_flux (((vector []){{taux.x,taux.y,taux.z},{{-1},{-1},{-1}}}));
     { 
#undef _OMPSTART
#define _OMPSTART double _maxres = maxres; 
#undef _OMPEND
#define _OMPEND OMP(omp critical) if (_maxres > maxres) maxres = _maxres; mpi_all_reduce_double (maxres, MPI_MAX); 
#line 113

if (!is_constant(rho)) {
#undef val_rho
#define val_rho(a,j,k,i) val(a,j,k,i)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,j,k,i)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,j,k,i)
#line 113
foreach (){

#line 113 "/home/popinet/basilisk-octree/src/viscosity.h"
 {

      double d = 0.;
      d += val(taux.z,0,0,1) - val(taux.z,0,0,0);

        d += val(taux.x,1,0,0) - val(taux.x,0,0,0);


        d += val(taux.y,0,1,0) - val(taux.y,0,0,0);

      val(res.z,0,0,0) = val(r.z,0,0,0) - ((coord){1.,1.,1.}).z*val(u.z,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.z,0,0,0)) > _maxres)
 _maxres = fabs (val(res.z,0,0,0));
    } } end_foreach(); }
if (is_constant(rho)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,j,k,i) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#line 113
foreach (){

#line 113 "/home/popinet/basilisk-octree/src/viscosity.h"
 {

      double d = 0.;
      d += val(taux.z,0,0,1) - val(taux.z,0,0,0);

        d += val(taux.x,1,0,0) - val(taux.x,0,0,0);


        d += val(taux.y,0,1,0) - val(taux.y,0,0,0);

      val(res.z,0,0,0) = val(r.z,0,0,0) - ((coord){1.,1.,1.}).z*val(u.z,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.z,0,0,0)) > _maxres)
 _maxres = fabs (val(res.z,0,0,0));
    } } end_foreach(); }
#undef _OMPSTART
#undef _OMPEND
#define _OMPSTART
#define _OMPEND
#line 127
 }
   delete (((scalar []){taux.x,taux.y,taux.z,{-1}})); }}
#line 157 "/home/popinet/basilisk-octree/src/viscosity.h"
  boundary (resl);
  return maxres;
}



mgstats viscosity (struct Viscosity p)
{
  vector u = p.u, r= new_vector("r");
   { foreach(){

#line 166 "/home/popinet/basilisk-octree/src/viscosity.h"

    {
#line 167

      val(r.x,0,0,0) = val(u.x,0,0,0);
#line 167

      val(r.y,0,0,0) = val(u.y,0,0,0);
#line 167

      val(r.z,0,0,0) = val(u.z,0,0,0);}; } end_foreach(); }

  vector mu = p.mu;
  scalar rho = p.rho;
  restriction (((scalar []){mu.x,mu.y,mu.z,rho,{-1}}));

  { mgstats _ret =  mg_solve ((scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}), (scalar *)((vector []){{r.x,r.y,r.z},{{-1},{-1},{-1}}}),
     residual_viscosity, relax_viscosity, &p); delete (((scalar []){r.x,r.y,r.z,{-1}}));  return _ret; }
 delete (((scalar []){r.x,r.y,r.z,{-1}})); }
#line 28 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
#line 37 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
scalar p= {0};
vector u= {{1},{2},{3}}, g= {{4},{5},{6}};
scalar pf= {7};
vector uf= {{8},{9},{10}};
#line 63 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
 vector mu = {{_NVARMAX + 0},{_NVARMAX + 1},{_NVARMAX + 2}}, a = {{_NVARMAX + 0},{_NVARMAX + 1},{_NVARMAX + 2}};
 vector alpha = {{_NVARMAX + 3},{_NVARMAX + 4},{_NVARMAX + 5}};
 scalar rho = {(_NVARMAX + 6)};
mgstats mgp, mgpf, mgu;
bool stokes = false;
#line 78 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
static void _set_boundary0 (void) { _attribute[p.i].boundary[right] = _boundary0; _attribute[p.i].boundary_homogeneous[right] = _boundary0_homogeneous; } 
static void _set_boundary1 (void) { _attribute[p.i].boundary[left] = _boundary1; _attribute[p.i].boundary_homogeneous[left] = _boundary1_homogeneous; } 





static void _set_boundary2 (void) { _attribute[p.i].boundary[top] = _boundary2; _attribute[p.i].boundary_homogeneous[top] = _boundary2_homogeneous; } 
static void _set_boundary3 (void) { _attribute[p.i].boundary[bottom] = _boundary3; _attribute[p.i].boundary_homogeneous[bottom] = _boundary3_homogeneous; } 


static void _set_boundary4 (void) { _attribute[p.i].boundary[front] = _boundary4; _attribute[p.i].boundary_homogeneous[front] = _boundary4_homogeneous; } 
static void _set_boundary5 (void) { _attribute[p.i].boundary[back] = _boundary5; _attribute[p.i].boundary_homogeneous[back] = _boundary5_homogeneous; } 






static int defaults_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults (const int i, const double t, Event * _ev) { trace ("defaults", "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 97); 
{




  CFL = 0.8;
   { foreach(){

#line 104 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
 {
    {
#line 105

      val(u.x,0,0,0) = val(g.x,0,0,0) = 0.;
#line 105

      val(u.y,0,0,0) = val(g.y,0,0,0) = 0.;
#line 105

      val(u.z,0,0,0) = val(g.z,0,0,0) = 0.;}
    val(p,0,0,0) = val(pf,0,0,0) = 0.;
  } } end_foreach(); }
   { foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 109
{

#line 109 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) = 0.; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 109
{

#line 109 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) = 0.; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 109
{

#line 109 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) = 0.; }  }}  end_foreach_face_generic()
#line 110
 end_foreach_face(); }




  if (!is_constant(alpha.x)) {
    vector alphav = alpha;
     { foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 117
{

#line 117 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

      val(alphav.x,0,0,0) = 1.; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 117
{

#line 117 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

      val(alphav.y,0,0,0) = 1.; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 117
{

#line 117 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

      val(alphav.z,0,0,0) = 1.; }  }}  end_foreach_face_generic()
#line 118
 end_foreach_face(); }
    boundary ((scalar *)((vector []){{alpha.x,alpha.y,alpha.z},{{-1},{-1},{-1}}}));
  }




  if (!is_constant(a.x)) {
    vector av = a;
     { foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 127
{

#line 127 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

      val(av.x,0,0,0) = 0.; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 127
{

#line 127 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

      val(av.y,0,0,0) = 0.; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 127
{

#line 127 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

      val(av.z,0,0,0) = 0.; }  }}  end_foreach_face_generic()
#line 128
 end_foreach_face(); }
    boundary ((scalar *)((vector []){{av.x,av.y,av.z},{{-1},{-1},{-1}}}));
  }

  boundary (((scalar []){p,pf,u.x,u.y,u.z,g.x,g.y,g.z,uf.x,uf.y,uf.z,{-1}}));






  _attribute[uf.x.i].refine = refine_face_solenoidal;

 end_trace("defaults", "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 141); } return 0; } 





static int init_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int init (const int i, const double t, Event * _ev) { trace ("init", "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 147); 
{
  boundary ((scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}));
  trash (((vector []){{uf.x,uf.y,uf.z},{{-1},{-1},{-1}}}));
   { 
if (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 151
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 151
{

#line 151 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(val(u.x,0,0,0) + val(u.x,-1,0,0))/2.; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 151
{

#line 151 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(val(u.y,0,0,0) + val(u.y,0,-1,0))/2.; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 151
{

#line 151 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*(val(u.z,0,0,0) + val(u.z,0,0,-1))/2.; }  }}  end_foreach_face_generic()
#line 152
 end_foreach_face(); }
if (is_constant(fm.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 151
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 151
{

#line 151 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(val(u.x,0,0,0) + val(u.x,-1,0,0))/2.; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 151
{

#line 151 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(val(u.y,0,0,0) + val(u.y,0,-1,0))/2.; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 151
{

#line 151 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*(val(u.z,0,0,0) + val(u.z,0,0,-1))/2.; }  }}  end_foreach_face_generic()
#line 152
 end_foreach_face(); } }
  boundary ((scalar *)((vector []){{uf.x,uf.y,uf.z},{{-1},{-1},{-1}}}));




  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  }

  event ("properties");
 end_trace("init", "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 164); } return 0; } 
#line 173 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
double dtmax;

static int set_dtmax_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int set_dtmax (const int i, const double t, Event * _ev) { trace ("set_dtmax", "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 175);  dtmax = DT; end_trace("set_dtmax", "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 175);  return 0; } 

static int stability_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int stability (const int i, const double t, Event * _ev) { trace ("stability", "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 177);  {
  dt = dtnext (t, timestep (uf, dtmax));
 end_trace("stability", "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 179); } return 0; } 







static int vof_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int vof (const int i, const double t, Event * _ev) { trace ("vof", "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 187); ; end_trace("vof", "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 187);  return 0; } 
static int tracer_advection_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int tracer_advection (const int i, const double t, Event * _ev) { trace ("tracer_advection", "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 188); ; end_trace("tracer_advection", "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 188);  return 0; } 






static int properties_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int properties (const int i, const double t, Event * _ev) { trace ("properties", "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 195);  {
  boundary (((scalar []){alpha.x,alpha.y,alpha.z,mu.x,mu.y,mu.z,rho,{-1}}));
 end_trace("properties", "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 197); } return 0; } 
#line 209 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
void prediction()
{
  vector du;
  {
#line 212
 {
    scalar s = new_scalar("s");
    du.x = s;
  }
#line 212
 {
    scalar s = new_scalar("s");
    du.y = s;
  }
#line 212
 {
    scalar s = new_scalar("s");
    du.z = s;
  }}

  if (_attribute[u.x.i].gradient)
     { foreach(){

#line 218 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

      {
#line 219

        val(du.x,0,0,0) = _attribute[u.x.i].gradient (val(u.x,-1,0,0), val(u.x,0,0,0), val(u.x,1,0,0))/Delta;
#line 219

        val(du.y,0,0,0) = _attribute[u.y.i].gradient (val(u.y,0,-1,0), val(u.y,0,0,0), val(u.y,0,1,0))/Delta;
#line 219

        val(du.z,0,0,0) = _attribute[u.z.i].gradient (val(u.z,0,0,-1), val(u.z,0,0,0), val(u.z,0,0,1))/Delta;}; } end_foreach(); }
  else
     { foreach(){

#line 222 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

      {
#line 223

        val(du.x,0,0,0) = (val(u.x,1,0,0) - val(u.x,-1,0,0))/(2.*Delta);
#line 223

        val(du.y,0,0,0) = (val(u.y,0,1,0) - val(u.y,0,-1,0))/(2.*Delta);
#line 223

        val(du.z,0,0,0) = (val(u.z,0,0,1) - val(u.z,0,0,-1))/(2.*Delta);}; } end_foreach(); }
  boundary ((scalar *)((vector []){{du.x,du.y,du.z},{{-1},{-1},{-1}}}));

  trash (((vector []){{uf.x,uf.y,uf.z},{{-1},{-1},{-1}}}));
   { 
if (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 228
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 228
{

#line 228 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. +
      s*min(1., 1. - s*un)*val(du.x,i,0,0)*Delta/2.;

      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);


      double fzz = val(u.z,i,0,0) < 0. ? val(u.z,i,0,1) - val(u.z,i,0,0) : val(u.z,i,0,0) - val(u.z,i,0,-1);
      val(uf.x,0,0,0) -= dt*val(u.z,i,0,0)*fzz/(2.*Delta);

    val(uf.x,0,0,0) *= val_fm_x(fm.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 228
{

#line 228 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. +
      s*min(1., 1. - s*un)*val(du.y,0,i,0)*Delta/2.;

      double fyy = val(u.z,0,i,0) < 0. ? val(u.y,0,i,1) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,0,i,-1);
      val(uf.y,0,0,0) -= dt*val(u.z,0,i,0)*fyy/(2.*Delta);


      double fzz = val(u.x,0,i,0) < 0. ? val(u.x,1,i,0) - val(u.x,0,i,0) : val(u.x,0,i,0) - val(u.x,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fzz/(2.*Delta);

    val(uf.y,0,0,0) *= val_fm_y(fm.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 228
{

#line 228 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.z,0,0,0) + val(u.z,0,0,-1))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.z,0,0,0) = val(u.z,0,0,i) + (val(g.z,0,0,0) + val(g.z,0,0,-1))*dt/4. +
      s*min(1., 1. - s*un)*val(du.z,0,0,i)*Delta/2.;

      double fyy = val(u.x,0,0,i) < 0. ? val(u.z,1,0,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,-1,0,i);
      val(uf.z,0,0,0) -= dt*val(u.x,0,0,i)*fyy/(2.*Delta);


      double fzz = val(u.y,0,0,i) < 0. ? val(u.y,0,1,i) - val(u.y,0,0,i) : val(u.y,0,0,i) - val(u.y,0,-1,i);
      val(uf.z,0,0,0) -= dt*val(u.y,0,0,i)*fzz/(2.*Delta);

    val(uf.z,0,0,0) *= val_fm_z(fm.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 242
 end_foreach_face(); }
if (is_constant(fm.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 228
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 228
{

#line 228 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. +
      s*min(1., 1. - s*un)*val(du.x,i,0,0)*Delta/2.;

      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);


      double fzz = val(u.z,i,0,0) < 0. ? val(u.z,i,0,1) - val(u.z,i,0,0) : val(u.z,i,0,0) - val(u.z,i,0,-1);
      val(uf.x,0,0,0) -= dt*val(u.z,i,0,0)*fzz/(2.*Delta);

    val(uf.x,0,0,0) *= val_fm_x(fm.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 228
{

#line 228 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. +
      s*min(1., 1. - s*un)*val(du.y,0,i,0)*Delta/2.;

      double fyy = val(u.z,0,i,0) < 0. ? val(u.y,0,i,1) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,0,i,-1);
      val(uf.y,0,0,0) -= dt*val(u.z,0,i,0)*fyy/(2.*Delta);


      double fzz = val(u.x,0,i,0) < 0. ? val(u.x,1,i,0) - val(u.x,0,i,0) : val(u.x,0,i,0) - val(u.x,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fzz/(2.*Delta);

    val(uf.y,0,0,0) *= val_fm_y(fm.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 228
{

#line 228 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.z,0,0,0) + val(u.z,0,0,-1))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.z,0,0,0) = val(u.z,0,0,i) + (val(g.z,0,0,0) + val(g.z,0,0,-1))*dt/4. +
      s*min(1., 1. - s*un)*val(du.z,0,0,i)*Delta/2.;

      double fyy = val(u.x,0,0,i) < 0. ? val(u.z,1,0,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,-1,0,i);
      val(uf.z,0,0,0) -= dt*val(u.x,0,0,i)*fyy/(2.*Delta);


      double fzz = val(u.y,0,0,i) < 0. ? val(u.y,0,1,i) - val(u.y,0,0,i) : val(u.y,0,0,i) - val(u.y,0,-1,i);
      val(uf.z,0,0,0) -= dt*val(u.y,0,0,i)*fzz/(2.*Delta);

    val(uf.z,0,0,0) *= val_fm_z(fm.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 242
 end_foreach_face(); } }
  boundary ((scalar *)((vector []){{uf.x,uf.y,uf.z},{{-1},{-1},{-1}}}));

  delete ((scalar *)((vector []){{du.x,du.y,du.z},{{-1},{-1},{-1}}}));
}
#line 257 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
static int advection_term_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int advection_term (const int i, const double t, Event * _ev) { trace ("advection_term", "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 257); 
{
  if (!stokes) {
    prediction();
    mgpf = project (uf, pf, alpha, dt/2.);
    advection ((struct Advection){(scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}), uf, dt, (scalar *)((vector []){{g.x,g.y,g.z},{{-1},{-1},{-1}}})});
  }
 end_trace("advection_term", "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 264); } return 0; } 







static void correction (double dt)
{
   { foreach(){

#line 274 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    {
#line 275

      val(u.x,0,0,0) += dt*val(g.x,0,0,0);
#line 275

      val(u.y,0,0,0) += dt*val(g.y,0,0,0);
#line 275

      val(u.z,0,0,0) += dt*val(g.z,0,0,0);}; } end_foreach(); }
  boundary ((scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}));
}
#line 287 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
static int viscous_term_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int viscous_term (const int i, const double t, Event * _ev) { trace ("viscous_term", "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 287); 
{
  if (constant(mu.x) != 0.) {
    correction (dt);
    mgu = viscosity ((struct Viscosity){u, mu, rho, dt});
    correction (-dt);
  }






  vector af = a;
  trash (((vector []){{uf.x,uf.y,uf.z},{af.x,af.y,af.z},{{-1},{-1},{-1}}}));
   { 
if (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 302
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 302
{

#line 302 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
 {
    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(val(u.x,0,0,0) + val(u.x,-1,0,0))/2.;
    if (!is_constant(af.x))
      val(af.x,0,0,0) = 0.;
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 302
{

#line 302 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
 {
    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(val(u.y,0,0,0) + val(u.y,0,-1,0))/2.;
    if (!is_constant(af.y))
      val(af.y,0,0,0) = 0.;
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 302
{

#line 302 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
 {
    val(uf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*(val(u.z,0,0,0) + val(u.z,0,0,-1))/2.;
    if (!is_constant(af.z))
      val(af.z,0,0,0) = 0.;
  } }  }}  end_foreach_face_generic()
#line 306
 end_foreach_face(); }
if (is_constant(fm.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 302
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 302
{

#line 302 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
 {
    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(val(u.x,0,0,0) + val(u.x,-1,0,0))/2.;
    if (!is_constant(af.x))
      val(af.x,0,0,0) = 0.;
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 302
{

#line 302 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
 {
    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(val(u.y,0,0,0) + val(u.y,0,-1,0))/2.;
    if (!is_constant(af.y))
      val(af.y,0,0,0) = 0.;
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 302
{

#line 302 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
 {
    val(uf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*(val(u.z,0,0,0) + val(u.z,0,0,-1))/2.;
    if (!is_constant(af.z))
      val(af.z,0,0,0) = 0.;
  } }  }}  end_foreach_face_generic()
#line 306
 end_foreach_face(); } }
 end_trace("viscous_term", "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 307); } return 0; } 
#line 322 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
static int acceleration_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int acceleration (const int i, const double t, Event * _ev) { trace ("acceleration", "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 322); 
{
  boundary_flux (((vector []){{a.x,a.y,a.z},{{-1},{-1},{-1}}}));
   { 
if (!is_constant(fm.x) && !is_constant(a.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#line 325
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 325
{

#line 325 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) += dt*val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 325
{

#line 325 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) += dt*val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 325
{

#line 325 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) += dt*val_fm_z(fm.z,0,0,0)*val_a_z(a.z,0,0,0); }  }}  end_foreach_face_generic()
#line 326
 end_foreach_face(); }
if (is_constant(fm.x) && !is_constant(a.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#line 325
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 325
{

#line 325 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) += dt*val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 325
{

#line 325 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) += dt*val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 325
{

#line 325 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) += dt*val_fm_z(fm.z,0,0,0)*val_a_z(a.z,0,0,0); }  }}  end_foreach_face_generic()
#line 326
 end_foreach_face(); }
if (!is_constant(fm.x) && is_constant(a.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#line 325
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 325
{

#line 325 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) += dt*val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 325
{

#line 325 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) += dt*val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 325
{

#line 325 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) += dt*val_fm_z(fm.z,0,0,0)*val_a_z(a.z,0,0,0); }  }}  end_foreach_face_generic()
#line 326
 end_foreach_face(); }
if (is_constant(fm.x) && is_constant(a.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#line 325
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 325
{

#line 325 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) += dt*val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 325
{

#line 325 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) += dt*val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 325
{

#line 325 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) += dt*val_fm_z(fm.z,0,0,0)*val_a_z(a.z,0,0,0); }  }}  end_foreach_face_generic()
#line 326
 end_foreach_face(); } }
  boundary ((scalar *)((vector []){{uf.x,uf.y,uf.z},{{-1},{-1},{-1}}}));
 end_trace("acceleration", "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 328); } return 0; } 
#line 337 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
static int projection_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int projection (const int i, const double t, Event * _ev) { trace ("projection", "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 337); 
{
  boundary (((scalar []){p,{-1}}));
  mgp = project (uf, p, alpha, dt);





  vector gf= new_face_vector("gf");
   { 
if (!is_constant(a.x) && !is_constant(alpha.x) && !is_constant(fm.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 347
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 348
 end_foreach_face(); }
if (is_constant(a.x) && !is_constant(alpha.x) && !is_constant(fm.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 347
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 348
 end_foreach_face(); }
if (!is_constant(a.x) && is_constant(alpha.x) && !is_constant(fm.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 347
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 348
 end_foreach_face(); }
if (is_constant(a.x) && is_constant(alpha.x) && !is_constant(fm.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 347
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 348
 end_foreach_face(); }
if (!is_constant(a.x) && !is_constant(alpha.x) && is_constant(fm.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 347
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 348
 end_foreach_face(); }
if (is_constant(a.x) && !is_constant(alpha.x) && is_constant(fm.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 347
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 348
 end_foreach_face(); }
if (!is_constant(a.x) && is_constant(alpha.x) && is_constant(fm.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 347
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 348
 end_foreach_face(); }
if (is_constant(a.x) && is_constant(alpha.x) && is_constant(fm.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 347
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 347
{

#line 347 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 348
 end_foreach_face(); } }
  boundary_flux (((vector []){{gf.x,gf.y,gf.z},{{-1},{-1},{-1}}}));





  trash (((vector []){{g.x,g.y,g.z},{{-1},{-1},{-1}}}));
   { foreach(){

#line 356 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

    {
#line 357

      val(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/2.;
#line 357

      val(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/2.;
#line 357

      val(g.z,0,0,0) = (val(gf.z,0,0,0) + val(gf.z,0,0,1))/2.;}; } end_foreach(); }
  boundary ((scalar *)((vector []){{g.x,g.y,g.z},{{-1},{-1},{-1}}}));




  correction (dt);
 delete (((scalar []){gf.x,gf.y,gf.z,{-1}}));  end_trace("projection", "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 365); } return 0; } 
#line 3 "atomisation.c"
#line 1 "vof.h"
#line 1 "/home/popinet/basilisk-octree/src/vof.h"
#line 15 "/home/popinet/basilisk-octree/src/vof.h"
#line 1 "fractions.h"
#line 1 "/home/popinet/basilisk-octree/src/fractions.h"
#line 11 "/home/popinet/basilisk-octree/src/fractions.h"
#line 1 "geometry.h"
#line 1 "/home/popinet/basilisk-octree/src/geometry.h"
#line 28 "/home/popinet/basilisk-octree/src/geometry.h"
double line_alpha (double c, coord n)
{
  double alpha, n1, n2;

  n1 = fabs (n.x); n2 = fabs (n.y);
  if (n1 > n2)
    swap (double, n1, n2);

  double v1 = n1/2.;
  if (c <= v1/n2)
    alpha = sqrt (2.*c*n1*n2);
  else if (c <= 1. - v1/n2)
    alpha = c*n2 + v1;
  else
    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;

  return alpha - (n.x + n.y)/2.;
}



double plane_alpha (double c, coord n)
{
  double alpha;
  coord n1;

  n1.x = fabs (n.x); n1.y = fabs (n.y); n1.z = fabs (n.z);

  double m1, m2, m3;
  m1 = min(n1.x, n1.y);
  m3 = max(n1.x, n1.y);
  m2 = n1.z;
  if (m2 < m1) {
    double tmp = m1;
    m1 = m2;
    m2 = tmp;
  }
  else if (m2 > m3) {
    double tmp = m3;
    m3 = m2;
    m2 = tmp;
  }
  double m12 = m1 + m2;
  double pr = max(6.*m1*m2*m3, 1e-50);
  double V1 = m1*m1*m1/pr;
  double V2 = V1 + (m2 - m1)/(2.*m3), V3;
  double mm;
  if (m3 < m12) {
    mm = m3;
    V3 = (m3*m3*(3.*m12 - m3) + m1*m1*(m1 - 3.*m3) + m2*m2*(m2 - 3.*m3))/pr;
  }
  else {
    mm = m12;
    V3 = mm/(2.*m3);
  }

  double ch = min(c, 1. - c);
  if (ch < V1)
    alpha = pow (pr*ch, 1./3.);
  else if (ch < V2)
    alpha = (m1 + sqrt(m1*m1 + 8.*m2*m3*(ch - V1)))/2.;
  else if (ch < V3) {
    double p12 = sqrt (2.*m1*m2);
    double q = 3.*(m12 - 2.*m3*ch)/(4.*p12);
    double teta = acos(clamp(q,-1.,1.))/3.;
    double cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + m12;
  }
  else if (m12 < m3)
    alpha = m3*ch + mm/2.;
  else {
    double p = m1*(m2 + m3) + m2*m3 - 1./4., p12 = sqrt(p);
    double q = 3.*m1*m2*m3*(1./2. - ch)/(2.*p*p12);
    double teta = acos(clamp(q,-1.,1.))/3.;
    double cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + 1./2.;
  }
  if (c > 1./2.) alpha = 1. - alpha;

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;
  if (n.z < 0.)
    alpha += n.z;

  return alpha - (n.x + n.y + n.z)/2.;;
}
#line 131 "/home/popinet/basilisk-octree/src/geometry.h"
double line_area (double nx, double ny, double alpha)
{
  double a, v, area;

  alpha += (nx + ny)/2.;
  if (nx < 0.) {
    alpha -= nx;
    nx = - nx;
  }
  if (ny < 0.) {
    alpha -= ny;
    ny = - ny;
  }

  if (alpha <= 0.)
    return 0.;

  if (alpha >= nx + ny)
    return 1.;

  if (nx < 1e-10)
    area = alpha/ny;
  else if (ny < 1e-10)
    area = alpha/nx;
  else {
    v = sq(alpha);

    a = alpha - nx;
    if (a > 0.)
      v -= a*a;

    a = alpha - ny;
    if (a > 0.)
      v -= a*a;

    area = v/(2.*nx*ny);
  }

  return clamp (area, 0., 1.);
}



double plane_volume (coord n, double alpha)
{
  double al = alpha + (n.x + n.y + n.z)/2. +
    max(0., -n.x) + max(0., -n.y) + max(0., -n.z);
  if (al <= 0.)
    return 0.;
  double tmp = fabs(n.x) + fabs(n.y) + fabs(n.z);
  if (al >= tmp)
    return 1.;
  if (tmp < 1e-10)
    return 0.;
  double n1 = fabs(n.x)/tmp;
  double n2 = fabs(n.y)/tmp;
  double n3 = fabs(n.z)/tmp;
  al = max(0., min(1., al/tmp));
  double al0 = min(al, 1. - al);
  double b1 = min(n1, n2);
  double b3 = max(n1, n2);
  double b2 = n3;
  if (b2 < b1) {
    tmp = b1;
    b1 = b2;
    b2 = tmp;
  }
  else if (b2 > b3) {
    tmp = b3;
    b3 = b2;
    b2 = tmp;
  }
  double b12 = b1 + b2;
  double bm = min(b12, b3);
  double pr = max(6.*b1*b2*b3, 1e-50);
  if (al0 < b1)
    tmp = al0*al0*al0/pr;
  else if (al0 < b2)
    tmp = 0.5*al0*(al0 - b1)/(b2*b3) + b1*b1*b1/pr;
  else if (al0 < bm)
    tmp = (al0*al0*(3.*b12 - al0) + b1*b1*(b1 - 3.*al0) +
    b2*b2*(b2 - 3.*al0))/pr;
  else if (b12 < b3)
    tmp = (al0 - 0.5*bm)/b3;
  else
    tmp = (al0*al0*(3. - 2.*al0) + b1*b1*(b1 - 3.*al0) +
    b2*b2*(b2 - 3.*al0) + b3*b3*(b3 - 3.*al0))/pr;

  double volume = al <= 0.5 ? tmp : 1. - tmp;
  return clamp (volume, 0., 1.);
}
#line 235 "/home/popinet/basilisk-octree/src/geometry.h"
double rectangle_fraction (coord n, double alpha, coord a, coord b)
{
  coord n1;
  {
#line 238
 {
    alpha -= n.x*(b.x + a.x)/2.;
    n1.x = n.x*(b.x - a.x);
  }
#line 238
 {
    alpha -= n.y*(b.y + a.y)/2.;
    n1.y = n.y*(b.y - a.y);
  }
#line 238
 {
    alpha -= n.z*(b.z + a.z)/2.;
    n1.z = n.z*(b.z - a.z);
  }}
  return plane_volume (n1, alpha);
}
#line 256 "/home/popinet/basilisk-octree/src/geometry.h"
int facets (double c, coord n, double alpha,
     coord p[2])
{
  assert (3 == 2);
  int i = 0;
  if (c > 0. && c < 1.) {
    for (double s = -0.5; s <= 0.5; s += 1.)
      {
#line 263

 if (fabs (n.y) > 1e-4 && i < 2) {
   double a = (alpha - s*n.x)/n.y;
   if (a >= -0.5 && a <= 0.5) {
     p[i].x = s;
     p[i++].y = a;
   }
 }
#line 263

 if (fabs (n.z) > 1e-4 && i < 2) {
   double a = (alpha - s*n.y)/n.z;
   if (a >= -0.5 && a <= 0.5) {
     p[i].y = s;
     p[i++].z = a;
   }
 }
#line 263

 if (fabs (n.x) > 1e-4 && i < 2) {
   double a = (alpha - s*n.z)/n.x;
   if (a >= -0.5 && a <= 0.5) {
     p[i].z = s;
     p[i++].x = a;
   }
 }}
  }
  return i;
}





double line_length_center (coord m, double alpha, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  if (n.x < 0.) {
    alpha -= n.x;
    n.x = - n.x;
  }
  if (n.y < 0.) {
    alpha -= n.y;
    n.y = - n.y;
  }

  p->x = p->y = 0.;

  if (alpha <= 0. || alpha >= n.x + n.y)
    return 0.;

  {
#line 298

    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = (m.y < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }
#line 298

    if (n.y < 1e-4) {
      p->y = 0.;
      p->x = (m.x < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }}

  if (alpha >= n.x) {
    p->x += 1.;
    p->y += (alpha - n.x)/n.y;
  }
  else
    p->x += alpha/n.x;

  double ax = p->x, ay = p->y;
  if (alpha >= n.y) {
    p->y += 1.;
    ay -= 1.;
    p->x += (alpha - n.y)/n.x;
    ax -= (alpha - n.y)/n.x;
  }
  else {
    p->y += alpha/n.y;
    ay -= alpha/n.y;
  }

  {
#line 324
 {
    p->x /= 2.;
    p->x = clamp (p->x, 0., 1.);
    if (m.x < 0.)
      p->x = 1. - p->x;
    p->x -= 0.5;
  }
#line 324
 {
    p->y /= 2.;
    p->y = clamp (p->y, 0., 1.);
    if (m.y < 0.)
      p->y = 1. - p->y;
    p->y -= 0.5;
  }}

  return sqrt (ax*ax + ay*ay);
}




double plane_area_center (coord m, double alpha, coord * p)
{
  if (fabs (m.x) < 1e-4) {
    coord n, q;
    n.x = m.y;
    n.y = m.z;
    double length = line_length_center (n, alpha, &q);
    p->x = 0.;
    p->y = q.x;
    p->z = q.y;
    return sq(length);
  }
  if (fabs (m.y) < 1e-4) {
    coord n, q;
    n.x = m.z;
    n.y = m.x;
    double length = line_length_center (n, alpha, &q);
    p->x = q.y;
    p->y = 0.;
    p->z = q.x;
    return sq(length);
  }
  if (fabs (m.z) < 1e-4) {
    double length = line_length_center (m, alpha, p);
    p->z = 0.;
    return sq(length);
  }

  alpha += (m.x + m.y + m.z)/2.;
  coord n = m;
  {
#line 368

    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
#line 368

    if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }
#line 368

    if (n.z < 0.) {
      alpha -= n.z;
      n.z = - n.z;
    }}

  double amax = n.x + n.y + n.z;
  if (alpha <= 0. || alpha >= amax) {
    p->x = p->y = p->z = 0.;
    return 0.;
  }

  double area = sq(alpha);
  p->x = p->y = p->z = area*alpha;

  {
#line 383
 {
    double b = alpha - n.x;
    if (b > 0.) {
      area -= b*b;
      p->x -= b*b*(2.*n.x + alpha);
      p->y -= b*b*b;
      p->z -= b*b*b;
    }
  }
#line 383
 {
    double b = alpha - n.y;
    if (b > 0.) {
      area -= b*b;
      p->y -= b*b*(2.*n.y + alpha);
      p->z -= b*b*b;
      p->x -= b*b*b;
    }
  }
#line 383
 {
    double b = alpha - n.z;
    if (b > 0.) {
      area -= b*b;
      p->z -= b*b*(2.*n.z + alpha);
      p->x -= b*b*b;
      p->y -= b*b*b;
    }
  }}

  amax = alpha - amax;
  {
#line 394
 {
    double b = amax + n.x;
    if (b > 0.) {
      area += b*b;
      p->y += b*b*(2.*n.y + alpha - n.z);
      p->z += b*b*(2.*n.z + alpha - n.y);
      p->x += b*b*b;
    }
  }
#line 394
 {
    double b = amax + n.y;
    if (b > 0.) {
      area += b*b;
      p->z += b*b*(2.*n.z + alpha - n.x);
      p->x += b*b*(2.*n.x + alpha - n.z);
      p->y += b*b*b;
    }
  }
#line 394
 {
    double b = amax + n.z;
    if (b > 0.) {
      area += b*b;
      p->x += b*b*(2.*n.x + alpha - n.y);
      p->y += b*b*(2.*n.y + alpha - n.x);
      p->z += b*b*b;
    }
  }}

  area *= 3.;
  {
#line 405
 {
    p->x /= area*n.x;
    p->x = clamp (p->x, 0., 1.);
    if (m.x < 0.) p->x = 1. - p->x;
    p->x -= 0.5;
  }
#line 405
 {
    p->y /= area*n.y;
    p->y = clamp (p->y, 0., 1.);
    if (m.y < 0.) p->y = 1. - p->y;
    p->y -= 0.5;
  }
#line 405
 {
    p->z /= area*n.z;
    p->z = clamp (p->z, 0., 1.);
    if (m.z < 0.) p->z = 1. - p->z;
    p->z -= 0.5;
  }}

  return area*sqrt (1./(sq(n.x)*sq(n.y)) +
      1./(sq(n.x)*sq(n.z)) +
      1./(sq(n.z)*sq(n.y)))/6.;
}
#line 12 "/home/popinet/basilisk-octree/src/fractions.h"
#line 20 "/home/popinet/basilisk-octree/src/fractions.h"
#line 1 "myc.h"
#line 1 "/home/popinet/basilisk-octree/src/myc.h"
#line 18 "/home/popinet/basilisk-octree/src/myc.h"
coord mycs (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 19 "/home/popinet/basilisk-octree/src/myc.h"

  double m1,m2,m[4][3],t0,t1,t2;
  int cn;



  m1 = val(c,-1,0,-1) + val(c,-1,0,1) + val(c,-1,-1,0) + val(c,-1,1,0) +
       val(c,-1,0,0);
  m2 = val(c,1,0,-1) + val(c,1,0,1) + val(c,1,-1,0) + val(c,1,1,0) +
       val(c,1,0,0);
  m[0][0] = m1 > m2 ? 1. : -1.;

  m1 = val(c,-1,-1,0)+ val(c,1,-1,0)+ val(c,0,-1,0);
  m2 = val(c,-1,1,0)+ val(c,1,1,0)+ val(c,0,1,0);
  m[0][1] = 0.5*(m1-m2);

  m1 = val(c,-1,0,-1)+ val(c,1,0,-1)+ val(c,0,0,-1);
  m2 = val(c,-1,0,1)+ val(c,1,0,1)+ val(c,0,0,1);
  m[0][2] = 0.5*(m1-m2);



  m1 = val(c,-1,-1,0) + val(c,-1,1,0) + val(c,-1,0,0);
  m2 = val(c,1,-1,0) + val(c,1,1,0) + val(c,1,0,0);
  m[1][0] = 0.5*(m1-m2);

  m1 = val(c,0,-1,-1) + val(c,0,-1,1) + val(c,1,-1,0) + val(c,-1,-1,0) +
       val(c,0,-1,0);
  m2 = val(c,0,1,-1) + val(c,0,1,1) + val(c,1,1,0) + val(c,-1,1,0) +
       val(c,0,1,0);
  m[1][1] = m1 > m2 ? 1. : -1.;

  m1 = val(c,0,-1,-1)+ val(c,0,0,-1)+ val(c,0,1,-1);
  m2 = val(c,0,-1,1)+ val(c,0,0,1)+ val(c,0,1,1);
  m[1][2] = 0.5*(m1-m2);




  m1 = val(c,-1,0,-1)+ val(c,-1,0,1)+ val(c,-1,0,0);
  m2 = val(c,1,0,-1)+ val(c,1,0,1)+ val(c,1,0,0);
  m[2][0] = 0.5*(m1-m2);

  m1 = val(c,0,-1,-1)+ val(c,0,-1,1)+ val(c,0,-1,0);
  m2 = val(c,0,1,-1)+ val(c,0,1,1)+ val(c,0,1,0);
  m[2][1] = 0.5*(m1-m2);

  m1 = val(c,-1,0,-1) + val(c,1,0,-1) + val(c,0,-1,-1) + val(c,0,1,-1) +
       val(c,0,0,-1);
  m2 = val(c,-1,0,1) + val(c,1,0,1) + val(c,0,-1,1) + val(c,0,1,1) +
       val(c,0,0,1);
  m[2][2] = m1 > m2 ? 1. : -1.;


  t0 = fabs(m[0][0]) + fabs(m[0][1]) + fabs(m[0][2]);
  m[0][0] /= t0;
  m[0][1] /= t0;
  m[0][2] /= t0;

  t0 = fabs(m[1][0]) + fabs(m[1][1]) + fabs(m[1][2]);
  m[1][0] /= t0;
  m[1][1] /= t0;
  m[1][2] /= t0;

  t0 = fabs(m[2][0]) + fabs(m[2][1]) + fabs(m[2][2]);
  m[2][0] /= t0;
  m[2][1] /= t0;
  m[2][2] /= t0;


  t0 = fabs(m[0][0]);
  t1 = fabs(m[1][1]);
  t2 = fabs(m[2][2]);

  cn = 0;
  if (t1 > t0) {
    t0 = t1;
    cn = 1;
  }
  if (t2 > t0)
    cn = 2;


  m1 = val(c,-1,-1,-1) + val(c,-1,1,-1) + val(c,-1,-1,1) + val(c,-1,1,1) +
       2.*(val(c,-1,-1,0) + val(c,-1,1,0) + val(c,-1,0,-1) + val(c,-1,0,1)) +
       4.*val(c,-1,0,0);
  m2 = val(c,1,-1,-1) + val(c,1,1,-1) + val(c,1,-1,1) + val(c,1,1,1) +
       2.*(val(c,1,-1,0) + val(c,1,1,0) + val(c,1,0,-1) + val(c,1,0,1)) +
       4.*val(c,1,0,0);
  m[3][0] = m1 - m2 + 1.e-30;

  m1 = val(c,-1,-1,-1) + val(c,-1,-1,1) + val(c,1,-1,-1) + val(c,1,-1,1) +
       2.*( val(c,-1,-1,0) + val(c,1,-1,0) + val(c,0,-1,-1) + val(c,0,-1,1)) +
       4.*val(c,0,-1,0);
  m2 = val(c,-1,1,-1) + val(c,-1,1,1) + val(c,1,1,-1) + val(c,1,1,1) +
       2.*(val(c,-1,1,0) + val(c,1,1,0) + val(c,0,1,-1) + val(c,0,1,1)) +
       4.*val(c,0,1,0);
  m[3][1] = m1 - m2 + 1.e-30;

  m1 = val(c,-1,-1,-1) + val(c,-1,1,-1) + val(c,1,-1,-1) + val(c,1,1,-1) +
       2.*(val(c,-1,0,-1) + val(c,1,0,-1) + val(c,0,-1,-1) + val(c,0,1,-1)) +
       4.*val(c,0,0,-1);
  m2 = val(c,-1,-1,1) + val(c,-1,1,1) + val(c,1,-1,1) + val(c,1,1,1) +
       2.*(val(c,-1,0,1) + val(c,1,0,1) + val(c,0,-1,1) + val(c,0,1,1)) +
       4.*val(c,0,0,1);
  m[3][2] = m1 - m2 + 1.e-30;


  t0 = fabs(m[3][0]) + fabs(m[3][1]) + fabs(m[3][2]);
  m[3][0] /= t0;
  m[3][1] /= t0;
  m[3][2] /= t0;


  t0 = fabs (m[3][0]);
  t1 = fabs (m[3][1]);
  t2 = fabs (m[3][2]);
  if (t1 > t0)
    t0 = t1;
  if (t2 > t0)
    t0 = t2;

  if (fabs(m[cn][cn]) > t0)
    cn = 3;


  coord mxyz;
  mxyz.x = m[cn][0];
  mxyz.y = m[cn][1];
  mxyz.z = m[cn][2];

  return mxyz;
}
#line 21 "/home/popinet/basilisk-octree/src/fractions.h"
#line 32 "/home/popinet/basilisk-octree/src/fractions.h"
void fraction_refine (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 33 "/home/popinet/basilisk-octree/src/fractions.h"






  double cc = val(c,0,0,0);
  if (cc <= 0. || cc >= 1.)
     { foreach_child()
      val(c,0,0,0) = cc; end_foreach_child(); }
  else {




    coord n = mycs (point, c);
    double alpha = plane_alpha (cc, n);






    coord a, b;
    {
#line 57
 {
      a.x = 0.; b.x = 0.5;
    }
#line 57
 {
      a.y = 0.; b.y = 0.5;
    }
#line 57
 {
      a.z = 0.; b.z = 0.5;
    }}

     { foreach_child() {
      coord nc;
      {
#line 63

 nc.x = child.x*n.x;
#line 63

 nc.y = child.y*n.y;
#line 63

 nc.z = child.z*n.z;}
      val(c,0,0,0) = rectangle_fraction (nc, alpha, a, b);
    } end_foreach_child(); }
  }
}











static void alpha_refine (Point point, scalar alpha)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 81 "/home/popinet/basilisk-octree/src/fractions.h"

  vector n = _attribute[alpha.i].n;
  double alphac = 2.*val(alpha,0,0,0);
  coord m;
  {
#line 85

    m.x = val(n.x,0,0,0);
#line 85

    m.y = val(n.y,0,0,0);
#line 85

    m.z = val(n.z,0,0,0);}
   { foreach_child() {
    val(alpha,0,0,0) = alphac;
    {
#line 89

      val(alpha,0,0,0) -= child.x*m.x/2.;
#line 89

      val(alpha,0,0,0) -= child.y*m.y/2.;
#line 89

      val(alpha,0,0,0) -= child.z*m.z/2.;}
  } end_foreach_child(); }
}
#line 116 "/home/popinet/basilisk-octree/src/fractions.h"
struct Fractions {
  scalar Phi;
  scalar c;
  vector s;
};

void fractions (struct Fractions a)
{
  scalar Phi = a.Phi;
  scalar c = a.c;
  vector s = (a.s).x.i ? (a.s) : new_face_vector("s");
#line 135 "/home/popinet/basilisk-octree/src/fractions.h"
  vector p= new_vector("p");
#line 147 "/home/popinet/basilisk-octree/src/fractions.h"
   { foreach_vertex(){

#line 147 "/home/popinet/basilisk-octree/src/fractions.h"
 {
#line 147
 if (neighbor(1,0,0).flags & vertex) {





    if (val(Phi,0,0,0)*val(Phi,1,0,0) < 0.) {






      val(p.x,0,0,0) = val(Phi,0,0,0)/(val(Phi,0,0,0) - val(Phi,1,0,0));
      if (val(Phi,0,0,0) < 0.)
 val(p.x,0,0,0) = 1. - val(p.x,0,0,0);
    }
#line 172 "/home/popinet/basilisk-octree/src/fractions.h"
    else
      val(p.x,0,0,0) = (val(Phi,0,0,0) > 0. || val(Phi,1,0,0) > 0.);
  }
#line 147
 if (neighbor(0,1,0).flags & vertex) {





    if (val(Phi,0,0,0)*val(Phi,0,1,0) < 0.) {






      val(p.y,0,0,0) = val(Phi,0,0,0)/(val(Phi,0,0,0) - val(Phi,0,1,0));
      if (val(Phi,0,0,0) < 0.)
 val(p.y,0,0,0) = 1. - val(p.y,0,0,0);
    }
#line 172 "/home/popinet/basilisk-octree/src/fractions.h"
    else
      val(p.y,0,0,0) = (val(Phi,0,0,0) > 0. || val(Phi,0,1,0) > 0.);
  }
#line 147
 if (neighbor(0,0,1).flags & vertex) {





    if (val(Phi,0,0,0)*val(Phi,0,0,1) < 0.) {






      val(p.z,0,0,0) = val(Phi,0,0,0)/(val(Phi,0,0,0) - val(Phi,0,0,1));
      if (val(Phi,0,0,0) < 0.)
 val(p.z,0,0,0) = 1. - val(p.z,0,0,0);
    }
#line 172 "/home/popinet/basilisk-octree/src/fractions.h"
    else
      val(p.z,0,0,0) = (val(Phi,0,0,0) > 0. || val(Phi,0,0,1) > 0.);
  }} } end_foreach_vertex(); }
#line 186 "/home/popinet/basilisk-octree/src/fractions.h"
  scalar s_x = s.x, s_y = s.y, s_z = s.z;
   { foreach_face_generic() { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 187
{

#line 187 "/home/popinet/basilisk-octree/src/fractions.h"






  {
#line 225 "/home/popinet/basilisk-octree/src/fractions.h"
    coord n;
    double nn = 0.;
    {
#line 227
 {
      n.x = val(p.y,0,0,0) - val(p.y,1,0,0);
      nn += fabs(n.x);
    }
#line 227
 {
      n.y = val(p.x,0,0,0) - val(p.x,0,1,0);
      nn += fabs(n.y);
    }}





    if (nn == 0.)
      val(s_z,0,0,0) = val(p.x,0,0,0);
    else {





      {
#line 244

 n.x /= nn;
#line 244

 n.y /= nn;}






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
#line 254

   if (val(p.x,0,i,0) > 0. && val(p.x,0,i,0) < 1.) {
     double a = sign(val(Phi,0,i,0))*(val(p.x,0,i,0) - 0.5);
     alpha += n.x*a + n.y*(i - 0.5);
     ni++;
   }
#line 254

   if (val(p.y,i,0,0) > 0. && val(p.y,i,0,0) < 1.) {
     double a = sign(val(Phi,i,0,0))*(val(p.y,i,0,0) - 0.5);
     alpha += n.y*a + n.x*(i - 0.5);
     ni++;
   }}






      val(s_z,0,0,0) = line_area (n.x, n.y, alpha/ni);
    }
  } }  }}  { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 187
{

#line 187 "/home/popinet/basilisk-octree/src/fractions.h"






  {
#line 225 "/home/popinet/basilisk-octree/src/fractions.h"
    coord n;
    double nn = 0.;
    {
#line 227
 {
      n.y = val(p.z,0,0,0) - val(p.z,0,1,0);
      nn += fabs(n.y);
    }
#line 227
 {
      n.z = val(p.y,0,0,0) - val(p.y,0,0,1);
      nn += fabs(n.z);
    }}





    if (nn == 0.)
      val(s_x,0,0,0) = val(p.y,0,0,0);
    else {





      {
#line 244

 n.y /= nn;
#line 244

 n.z /= nn;}






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
#line 254

   if (val(p.y,0,0,i) > 0. && val(p.y,0,0,i) < 1.) {
     double a = sign(val(Phi,0,0,i))*(val(p.y,0,0,i) - 0.5);
     alpha += n.y*a + n.z*(i - 0.5);
     ni++;
   }
#line 254

   if (val(p.z,0,i,0) > 0. && val(p.z,0,i,0) < 1.) {
     double a = sign(val(Phi,0,i,0))*(val(p.z,0,i,0) - 0.5);
     alpha += n.z*a + n.y*(i - 0.5);
     ni++;
   }}






      val(s_x,0,0,0) = line_area (n.y, n.z, alpha/ni);
    }
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 187
{

#line 187 "/home/popinet/basilisk-octree/src/fractions.h"






  {
#line 225 "/home/popinet/basilisk-octree/src/fractions.h"
    coord n;
    double nn = 0.;
    {
#line 227
 {
      n.z = val(p.x,0,0,0) - val(p.x,0,0,1);
      nn += fabs(n.z);
    }
#line 227
 {
      n.x = val(p.z,0,0,0) - val(p.z,1,0,0);
      nn += fabs(n.x);
    }}





    if (nn == 0.)
      val(s_y,0,0,0) = val(p.z,0,0,0);
    else {





      {
#line 244

 n.z /= nn;
#line 244

 n.x /= nn;}






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
#line 254

   if (val(p.z,i,0,0) > 0. && val(p.z,i,0,0) < 1.) {
     double a = sign(val(Phi,i,0,0))*(val(p.z,i,0,0) - 0.5);
     alpha += n.z*a + n.x*(i - 0.5);
     ni++;
   }
#line 254

   if (val(p.x,0,0,i) > 0. && val(p.x,0,0,i) < 1.) {
     double a = sign(val(Phi,0,0,i))*(val(p.x,0,0,i) - 0.5);
     alpha += n.x*a + n.z*(i - 0.5);
     ni++;
   }}






      val(s_y,0,0,0) = line_area (n.z, n.x, alpha/ni);
    }
  } }  }}  end_foreach_face_generic()
#line 268
 end_foreach_face(); }







  boundary_flux (((vector []){{s.x,s.y,s.z},{{-1},{-1},{-1}}}));
   { foreach(){

#line 277 "/home/popinet/basilisk-octree/src/fractions.h"
 {




    coord n;
    double nn = 0.;
    {
#line 284
 {
      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);
      nn += fabs(n.x);
    }
#line 284
 {
      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);
      nn += fabs(n.y);
    }
#line 284
 {
      n.z = val(s.z,0,0,0) - val(s.z,0,0,1);
      nn += fabs(n.z);
    }}
    if (nn == 0.)
      val(c,0,0,0) = val(s.x,0,0,0);
    else {
      {
#line 291

 n.x /= nn;
#line 291

 n.y /= nn;
#line 291

 n.z /= nn;}






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)
   {
#line 302

     if (val(p.x,0,i,j) > 0. && val(p.x,0,i,j) < 1.) {
       double a = sign(val(Phi,0,i,j))*(val(p.x,0,i,j) - 0.5);
       alpha += n.x*a + n.y*(i - 0.5) + n.z*(j - 0.5);
       ni++;
     }
#line 302

     if (val(p.y,j,0,i) > 0. && val(p.y,j,0,i) < 1.) {
       double a = sign(val(Phi,j,0,i))*(val(p.y,j,0,i) - 0.5);
       alpha += n.y*a + n.z*(i - 0.5) + n.x*(j - 0.5);
       ni++;
     }
#line 302

     if (val(p.z,i,j,0) > 0. && val(p.z,i,j,0) < 1.) {
       double a = sign(val(Phi,i,j,0))*(val(p.z,i,j,0) - 0.5);
       alpha += n.z*a + n.x*(i - 0.5) + n.y*(j - 0.5);
       ni++;
     }}




      val(c,0,0,0) = ni ? plane_volume (n, alpha/ni) : val(s.x,0,0,0);
    }
  } } end_foreach(); }





  boundary (((scalar []){c,{-1}}));
 delete (((scalar []){p.x,p.y,p.z,{-1}}));  { if (!(a.s).x.i) delete (((scalar []){s.x,s.y,s.z,{-1}})); } }
#line 335 "/home/popinet/basilisk-octree/src/fractions.h"
coord youngs_normal (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 336 "/home/popinet/basilisk-octree/src/fractions.h"

  coord n;
  double nn = 0.;
  assert (3 == 2);
  {
#line 340
 {
    n.x = (val(c,-1,1,0) + 2.*val(c,-1,0,0) + val(c,-1,-1,0) -
    val(c,+1,1,0) - 2.*val(c,+1,0,0) - val(c,+1,-1,0));
    nn += fabs(n.x);
  }
#line 340
 {
    n.y = (val(c,0,-1,1) + 2.*val(c,0,-1,0) + val(c,0,-1,-1) -
    val(c,0,+1,1) - 2.*val(c,0,+1,0) - val(c,0,+1,-1));
    nn += fabs(n.y);
  }
#line 340
 {
    n.z = (val(c,1,0,-1) + 2.*val(c,0,0,-1) + val(c,-1,0,-1) -
    val(c,1,0,+1) - 2.*val(c,0,0,+1) - val(c,-1,0,+1));
    nn += fabs(n.z);
  }}

  if (nn > 0.)
    {
#line 347

      n.x /= nn;
#line 347

      n.y /= nn;
#line 347

      n.z /= nn;}
  else
    n.x = 1.;
  return n;
}
#line 361 "/home/popinet/basilisk-octree/src/fractions.h"
void reconstruction (const scalar c, vector n, scalar alpha)
{
   { foreach(){

#line 363 "/home/popinet/basilisk-octree/src/fractions.h"
 {





    if (val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) {
      val(alpha,0,0,0) = 0.;
      {
#line 371

 val(n.x,0,0,0) = 0.;
#line 371

 val(n.y,0,0,0) = 0.;
#line 371

 val(n.z,0,0,0) = 0.;}
    }
    else {






      coord m = mycs (point, c);

      {
#line 383

 val(n.x,0,0,0) = m.x;
#line 383

 val(n.y,0,0,0) = m.y;
#line 383

 val(n.z,0,0,0) = m.z;}
      val(alpha,0,0,0) = plane_alpha (val(c,0,0,0), m);
    }
  } } end_foreach(); }
#line 396 "/home/popinet/basilisk-octree/src/fractions.h"
  {
#line 396

    _attribute[n.x.i].refine = _attribute[n.x.i].prolongation = refine_injection;
#line 396

    _attribute[n.y.i].refine = _attribute[n.y.i].prolongation = refine_injection;
#line 396

    _attribute[n.z.i].refine = _attribute[n.z.i].prolongation = refine_injection;}




  _attribute[alpha.i].n = n;
  _attribute[alpha.i].refine = _attribute[alpha.i].prolongation = alpha_refine;







  boundary (((scalar []){n.x,n.y,n.z,alpha,{-1}}));
}
#line 432 "/home/popinet/basilisk-octree/src/fractions.h"
struct OutputFacets {
  scalar c;
  FILE * fp;
  vector s;
};

void output_facets (struct OutputFacets p)
{
  assert (3 == 2);

  scalar c = p.c;
  vector s = p.s;
  if (!p.fp) p.fp = qstdout();

   { foreach(){

#line 446 "/home/popinet/basilisk-octree/src/fractions.h"

    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n;
      if (!s.x.i)

 n = mycs (point, c);
      else {

 double nn = 0.;
 {
#line 455
 {
   n.x = val(s.x,0,0,0) - val(s.x,1,0,0);
   nn += fabs(n.x);
 }
#line 455
 {
   n.y = val(s.y,0,0,0) - val(s.y,0,1,0);
   nn += fabs(n.y);
 }
#line 455
 {
   n.z = val(s.z,0,0,0) - val(s.z,0,0,1);
   nn += fabs(n.z);
 }}
 assert (nn > 0.);
 {
#line 460

   n.x /= nn;
#line 460

   n.y /= nn;
#line 460

   n.z /= nn;}
      }
      double alpha = plane_alpha (val(c,0,0,0), n);
      coord segment[2];
      if (facets (val(c,0,0,0), n, alpha, segment) == 2)
 fprintf (p.fp, "%g %g\n%g %g\n\n",
   x + segment[0].x*Delta, y + segment[0].y*Delta,
   x + segment[1].x*Delta, y + segment[1].y*Delta);
    } } end_foreach(); }

  fflush (p.fp);
}
#line 16 "/home/popinet/basilisk-octree/src/vof.h"
#line 24 "/home/popinet/basilisk-octree/src/vof.h"
extern scalar * interfaces;
extern vector uf;
extern double dt;





static int defaults_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults_0 (const int i, const double t, Event * _ev) { trace ("defaults_0", "/home/popinet/basilisk-octree/src/vof.h", 32); 
{

  if (interfaces) for (scalar c = *interfaces, *_i79 = interfaces; ((scalar *)&c)->i >= 0; c = *++_i79)
    _attribute[c.i].refine = _attribute[c.i].prolongation = fraction_refine;

 end_trace("defaults_0", "/home/popinet/basilisk-octree/src/vof.h", 38); } return 0; } 





static int stability_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int stability_0 (const int i, const double t, Event * _ev) { trace ("stability_0", "/home/popinet/basilisk-octree/src/vof.h", 44);  {
  if (CFL > 0.5)
    CFL = 0.5;
 end_trace("stability_0", "/home/popinet/basilisk-octree/src/vof.h", 47); } return 0; } 
#line 61 "/home/popinet/basilisk-octree/src/vof.h"

#line 61

static void sweep_x (scalar c, scalar cc)
{
  vector n= new_vector("n");
  scalar alpha= new_scalar("alpha"), flux= new_scalar("flux");
  double cfl = 0.;






  reconstruction (c, n, alpha);

   { 
#undef _OMPSTART
#define _OMPSTART double _cfl = cfl; 
#undef _OMPEND
#define _OMPEND OMP(omp critical) if (_cfl > cfl) cfl = _cfl; mpi_all_reduce_double (cfl, MPI_MAX); 
#line 75

if (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 75
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 75
{

#line 75 "/home/popinet/basilisk-octree/src/vof.h"
 {






    double un = val(uf.x,0,0,0)*dt/(Delta*val_fm_x(fm.x,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0);
#line 102 "/home/popinet/basilisk-octree/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), val(n.z,i,0,0)}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); }
if (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 75
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 75
{

#line 75 "/home/popinet/basilisk-octree/src/vof.h"
 {






    double un = val(uf.x,0,0,0)*dt/(Delta*val_fm_x(fm.x,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0);
#line 102 "/home/popinet/basilisk-octree/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), val(n.z,i,0,0)}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); }
if (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 75
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 75
{

#line 75 "/home/popinet/basilisk-octree/src/vof.h"
 {






    double un = val(uf.x,0,0,0)*dt/(Delta*val_fm_x(fm.x,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0);
#line 102 "/home/popinet/basilisk-octree/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), val(n.z,i,0,0)}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); }
if (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 75
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 75
{

#line 75 "/home/popinet/basilisk-octree/src/vof.h"
 {






    double un = val(uf.x,0,0,0)*dt/(Delta*val_fm_x(fm.x,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_x(fm.x,0,0,0)*s/val_cm(cm,0,0,0);
#line 102 "/home/popinet/basilisk-octree/src/vof.h"
    double cf = (val(c,i,0,0) <= 0. || val(c,i,0,0) >= 1.) ? val(c,i,0,0) :
      rectangle_fraction ((coord){-s*val(n.x,i,0,0), val(n.y,i,0,0), val(n.z,i,0,0)}, val(alpha,i,0,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.x,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); }
#undef _OMPSTART
#undef _OMPEND
#define _OMPSTART
#define _OMPEND
#line 113
 }
#line 122 "/home/popinet/basilisk-octree/src/vof.h"
  for (int l = depth() - 1; l >= 0; l--)
     { foreach_halo (prolongation, l){

#line 123 "/home/popinet/basilisk-octree/src/vof.h"
 {
#line 135 "/home/popinet/basilisk-octree/src/vof.h"
      if ((!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0))
 val(flux,0,0,0) = (fine(flux,0,0,0) + fine(flux,0,1,0) +
    fine(flux,0,0,1) + fine(flux,0,1,1))/4.;
      if ((!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0))
 val(flux,1,0,0) = (fine(flux,2,0,0) + fine(flux,2,1,0) +
     fine(flux,2,0,1) + fine(flux,2,1,1))/4.;

    } } end_foreach_halo(); }





  if (cfl > 0.5 + 1e-6)
    fprintf (ferr,
      "WARNING: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n",
      cfl - 0.5), fflush (ferr);
#line 164 "/home/popinet/basilisk-octree/src/vof.h"
   { 
if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 164
foreach(){

#line 164 "/home/popinet/basilisk-octree/src/vof.h"

    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,1,0,0) + val(cc,0,0,0)*(val(uf.x,1,0,0) - val(uf.x,0,0,0)))/(val_cm(cm,0,0,0)*Delta); } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 164
foreach(){

#line 164 "/home/popinet/basilisk-octree/src/vof.h"

    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,1,0,0) + val(cc,0,0,0)*(val(uf.x,1,0,0) - val(uf.x,0,0,0)))/(val_cm(cm,0,0,0)*Delta); } end_foreach(); } }
  boundary (((scalar []){c,{-1}}));
 delete (((scalar []){flux,alpha,n.x,n.y,n.z,{-1}})); }
#line 61

static void sweep_y (scalar c, scalar cc)
{
  vector n= new_vector("n");
  scalar alpha= new_scalar("alpha"), flux= new_scalar("flux");
  double cfl = 0.;






  reconstruction (c, n, alpha);

   { 
#undef _OMPSTART
#define _OMPSTART double _cfl = cfl; 
#undef _OMPEND
#define _OMPEND OMP(omp critical) if (_cfl > cfl) cfl = _cfl; mpi_all_reduce_double (cfl, MPI_MAX); 
#line 75

if (!is_constant(fm.y) && !is_constant(cm)) {
#undef val_fm_y
#define val_fm_y(a,k,i,j) val(a,k,i,j)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,k,i,j)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,k,i,j)
#undef val_fm_z
#define val_fm_z(a,k,i,j) val(a,k,i,j)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,k,i,j)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,k,i,j)
#undef val_fm_x
#define val_fm_x(a,k,i,j) val(a,k,i,j)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,k,i,j)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,k,i,j)
#undef val_cm
#define val_cm(a,k,i,j) val(a,k,i,j)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,k,i,j)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,k,i,j)
#line 75
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_y()) {
#line 75
{

#line 75 "/home/popinet/basilisk-octree/src/vof.h"
 {






    double un = val(uf.y,0,0,0)*dt/(Delta*val_fm_y(fm.y,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0);
#line 102 "/home/popinet/basilisk-octree/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.z,0,i,0), val(n.x,0,i,0)}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); }
if (is_constant(fm.y) && !is_constant(cm)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.y.i -_NVARMAX], _constant[fm.z.i - _NVARMAX], _constant[fm.x.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_y
#define val_fm_y(a,k,i,j) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,k,i,j) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_fm_x
#define val_fm_x(a,k,i,j) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_cm
#define val_cm(a,k,i,j) val(a,k,i,j)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,k,i,j)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,k,i,j)
#line 75
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_y()) {
#line 75
{

#line 75 "/home/popinet/basilisk-octree/src/vof.h"
 {






    double un = val(uf.y,0,0,0)*dt/(Delta*val_fm_y(fm.y,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0);
#line 102 "/home/popinet/basilisk-octree/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.z,0,i,0), val(n.x,0,i,0)}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); }
if (!is_constant(fm.y) && is_constant(cm)) {
#undef val_fm_y
#define val_fm_y(a,k,i,j) val(a,k,i,j)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,k,i,j)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,k,i,j)
#undef val_fm_z
#define val_fm_z(a,k,i,j) val(a,k,i,j)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,k,i,j)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,k,i,j)
#undef val_fm_x
#define val_fm_x(a,k,i,j) val(a,k,i,j)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,k,i,j)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,k,i,j)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,k,i,j) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 75
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_y()) {
#line 75
{

#line 75 "/home/popinet/basilisk-octree/src/vof.h"
 {






    double un = val(uf.y,0,0,0)*dt/(Delta*val_fm_y(fm.y,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0);
#line 102 "/home/popinet/basilisk-octree/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.z,0,i,0), val(n.x,0,i,0)}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); }
if (is_constant(fm.y) && is_constant(cm)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.y.i -_NVARMAX], _constant[fm.z.i - _NVARMAX], _constant[fm.x.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_y
#define val_fm_y(a,k,i,j) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,k,i,j) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_fm_x
#define val_fm_x(a,k,i,j) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,k,i,j) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 75
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_y()) {
#line 75
{

#line 75 "/home/popinet/basilisk-octree/src/vof.h"
 {






    double un = val(uf.y,0,0,0)*dt/(Delta*val_fm_y(fm.y,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_y(fm.y,0,0,0)*s/val_cm(cm,0,0,0);
#line 102 "/home/popinet/basilisk-octree/src/vof.h"
    double cf = (val(c,0,i,0) <= 0. || val(c,0,i,0) >= 1.) ? val(c,0,i,0) :
      rectangle_fraction ((coord){-s*val(n.y,0,i,0), val(n.z,0,i,0), val(n.x,0,i,0)}, val(alpha,0,i,0),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.y,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); }
#undef _OMPSTART
#undef _OMPEND
#define _OMPSTART
#define _OMPEND
#line 113
 }
#line 122 "/home/popinet/basilisk-octree/src/vof.h"
  for (int l = depth() - 1; l >= 0; l--)
     { foreach_halo (prolongation, l){

#line 123 "/home/popinet/basilisk-octree/src/vof.h"
 {
#line 135 "/home/popinet/basilisk-octree/src/vof.h"
      if ((!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0))
 val(flux,0,0,0) = (fine(flux,0,0,0) + fine(flux,0,0,1) +
    fine(flux,1,0,0) + fine(flux,1,0,1))/4.;
      if ((!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0))
 val(flux,0,1,0) = (fine(flux,0,2,0) + fine(flux,0,2,1) +
     fine(flux,1,2,0) + fine(flux,1,2,1))/4.;

    } } end_foreach_halo(); }





  if (cfl > 0.5 + 1e-6)
    fprintf (ferr,
      "WARNING: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n",
      cfl - 0.5), fflush (ferr);
#line 164 "/home/popinet/basilisk-octree/src/vof.h"
   { 
if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,k,i,j) val(a,k,i,j)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,k,i,j)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,k,i,j)
#line 164
foreach(){

#line 164 "/home/popinet/basilisk-octree/src/vof.h"

    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,0,1,0) + val(cc,0,0,0)*(val(uf.y,0,1,0) - val(uf.y,0,0,0)))/(val_cm(cm,0,0,0)*Delta); } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,k,i,j) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 164
foreach(){

#line 164 "/home/popinet/basilisk-octree/src/vof.h"

    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,0,1,0) + val(cc,0,0,0)*(val(uf.y,0,1,0) - val(uf.y,0,0,0)))/(val_cm(cm,0,0,0)*Delta); } end_foreach(); } }
  boundary (((scalar []){c,{-1}}));
 delete (((scalar []){flux,alpha,n.x,n.y,n.z,{-1}})); }
#line 61

static void sweep_z (scalar c, scalar cc)
{
  vector n= new_vector("n");
  scalar alpha= new_scalar("alpha"), flux= new_scalar("flux");
  double cfl = 0.;






  reconstruction (c, n, alpha);

   { 
#undef _OMPSTART
#define _OMPSTART double _cfl = cfl; 
#undef _OMPEND
#define _OMPEND OMP(omp critical) if (_cfl > cfl) cfl = _cfl; mpi_all_reduce_double (cfl, MPI_MAX); 
#line 75

if (!is_constant(fm.z) && !is_constant(cm)) {
#undef val_fm_z
#define val_fm_z(a,j,k,i) val(a,j,k,i)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,j,k,i)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,j,k,i)
#undef val_fm_x
#define val_fm_x(a,j,k,i) val(a,j,k,i)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,j,k,i)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,j,k,i)
#undef val_fm_y
#define val_fm_y(a,j,k,i) val(a,j,k,i)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,j,k,i)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,j,k,i)
#undef val_cm
#define val_cm(a,j,k,i) val(a,j,k,i)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,j,k,i)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,j,k,i)
#line 75
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_z()) {
#line 75
{

#line 75 "/home/popinet/basilisk-octree/src/vof.h"
 {






    double un = val(uf.z,0,0,0)*dt/(Delta*val_fm_z(fm.z,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_z(fm.z,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_z(fm.z,0,0,0)*s/val_cm(cm,0,0,0);
#line 102 "/home/popinet/basilisk-octree/src/vof.h"
    double cf = (val(c,0,0,i) <= 0. || val(c,0,0,i) >= 1.) ? val(c,0,0,i) :
      rectangle_fraction ((coord){-s*val(n.z,0,0,i), val(n.x,0,0,i), val(n.y,0,0,i)}, val(alpha,0,0,i),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); }
if (is_constant(fm.z) && !is_constant(cm)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.z.i -_NVARMAX], _constant[fm.x.i - _NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_z
#define val_fm_z(a,j,k,i) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_fm_x
#define val_fm_x(a,j,k,i) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,j,k,i) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_cm
#define val_cm(a,j,k,i) val(a,j,k,i)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,j,k,i)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,j,k,i)
#line 75
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_z()) {
#line 75
{

#line 75 "/home/popinet/basilisk-octree/src/vof.h"
 {






    double un = val(uf.z,0,0,0)*dt/(Delta*val_fm_z(fm.z,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_z(fm.z,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_z(fm.z,0,0,0)*s/val_cm(cm,0,0,0);
#line 102 "/home/popinet/basilisk-octree/src/vof.h"
    double cf = (val(c,0,0,i) <= 0. || val(c,0,0,i) >= 1.) ? val(c,0,0,i) :
      rectangle_fraction ((coord){-s*val(n.z,0,0,i), val(n.x,0,0,i), val(n.y,0,0,i)}, val(alpha,0,0,i),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); }
if (!is_constant(fm.z) && is_constant(cm)) {
#undef val_fm_z
#define val_fm_z(a,j,k,i) val(a,j,k,i)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,j,k,i)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,j,k,i)
#undef val_fm_x
#define val_fm_x(a,j,k,i) val(a,j,k,i)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,j,k,i)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,j,k,i)
#undef val_fm_y
#define val_fm_y(a,j,k,i) val(a,j,k,i)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,j,k,i)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,j,k,i)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,j,k,i) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 75
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_z()) {
#line 75
{

#line 75 "/home/popinet/basilisk-octree/src/vof.h"
 {






    double un = val(uf.z,0,0,0)*dt/(Delta*val_fm_z(fm.z,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_z(fm.z,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_z(fm.z,0,0,0)*s/val_cm(cm,0,0,0);
#line 102 "/home/popinet/basilisk-octree/src/vof.h"
    double cf = (val(c,0,0,i) <= 0. || val(c,0,0,i) >= 1.) ? val(c,0,0,i) :
      rectangle_fraction ((coord){-s*val(n.z,0,0,i), val(n.x,0,0,i), val(n.y,0,0,i)}, val(alpha,0,0,i),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); }
if (is_constant(fm.z) && is_constant(cm)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.z.i -_NVARMAX], _constant[fm.x.i - _NVARMAX], _constant[fm.y.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_z
#define val_fm_z(a,j,k,i) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_fm_x
#define val_fm_x(a,j,k,i) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,j,k,i) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,j,k,i) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 75
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_z()) {
#line 75
{

#line 75 "/home/popinet/basilisk-octree/src/vof.h"
 {






    double un = val(uf.z,0,0,0)*dt/(Delta*val_fm_z(fm.z,0,0,0)), s = sign(un);
    int i = -(s + 1.)/2.;




    if (un*val_fm_z(fm.z,0,0,0)*s/val_cm(cm,0,0,0) > _cfl)
      _cfl = un*val_fm_z(fm.z,0,0,0)*s/val_cm(cm,0,0,0);
#line 102 "/home/popinet/basilisk-octree/src/vof.h"
    double cf = (val(c,0,0,i) <= 0. || val(c,0,0,i) >= 1.) ? val(c,0,0,i) :
      rectangle_fraction ((coord){-s*val(n.z,0,0,i), val(n.x,0,0,i), val(n.y,0,0,i)}, val(alpha,0,0,i),
     (coord){-0.5, -0.5, -0.5},
     (coord){s*un - 0.5, 0.5, 0.5});





    val(flux,0,0,0) = cf*val(uf.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); }
#undef _OMPSTART
#undef _OMPEND
#define _OMPSTART
#define _OMPEND
#line 113
 }
#line 122 "/home/popinet/basilisk-octree/src/vof.h"
  for (int l = depth() - 1; l >= 0; l--)
     { foreach_halo (prolongation, l){

#line 123 "/home/popinet/basilisk-octree/src/vof.h"
 {
#line 135 "/home/popinet/basilisk-octree/src/vof.h"
      if ((!is_leaf (neighbor(0,0,-1)) && neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0))
 val(flux,0,0,0) = (fine(flux,0,0,0) + fine(flux,1,0,0) +
    fine(flux,0,1,0) + fine(flux,1,1,0))/4.;
      if ((!is_leaf (neighbor(0,0,1)) && neighbor(0,0,1).neighbors && neighbor(0,0,1).pid >= 0))
 val(flux,0,0,1) = (fine(flux,0,0,2) + fine(flux,1,0,2) +
     fine(flux,0,1,2) + fine(flux,1,1,2))/4.;

    } } end_foreach_halo(); }





  if (cfl > 0.5 + 1e-6)
    fprintf (ferr,
      "WARNING: CFL must be <= 0.5 for VOF (cfl - 0.5 = %g)\n",
      cfl - 0.5), fflush (ferr);
#line 164 "/home/popinet/basilisk-octree/src/vof.h"
   { 
if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,j,k,i) val(a,j,k,i)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,j,k,i)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,j,k,i)
#line 164
foreach(){

#line 164 "/home/popinet/basilisk-octree/src/vof.h"

    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,0,0,1) + val(cc,0,0,0)*(val(uf.z,0,0,1) - val(uf.z,0,0,0)))/(val_cm(cm,0,0,0)*Delta); } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,j,k,i) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 164
foreach(){

#line 164 "/home/popinet/basilisk-octree/src/vof.h"

    val(c,0,0,0) += dt*(val(flux,0,0,0) - val(flux,0,0,1) + val(cc,0,0,0)*(val(uf.z,0,0,1) - val(uf.z,0,0,0)))/(val_cm(cm,0,0,0)*Delta); } end_foreach(); } }
  boundary (((scalar []){c,{-1}}));
 delete (((scalar []){flux,alpha,n.x,n.y,n.z,{-1}})); }






static int vof_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int vof_0 (const int i, const double t, Event * _ev) { trace ("vof_0", "/home/popinet/basilisk-octree/src/vof.h", 174); 
{
  if (interfaces) for (scalar c = *interfaces, *_i80 = interfaces; ((scalar *)&c)->i >= 0; c = *++_i80) {
#line 186 "/home/popinet/basilisk-octree/src/vof.h"
    scalar cc= new_scalar("cc");
     { foreach(){

#line 187 "/home/popinet/basilisk-octree/src/vof.h"

      val(cc,0,0,0) = (val(c,0,0,0) > 0.5); } end_foreach(); }






    void (* sweep[3]) (scalar, scalar);
    int d = 0;
    {
#line 197

      sweep[d++] = sweep_x;
#line 197

      sweep[d++] = sweep_y;
#line 197

      sweep[d++] = sweep_z;}
    boundary (((scalar []){c,{-1}}));
    for (d = 0; d < 3; d++)
      sweep[(i + d) % 3] (c, cc);
   delete (((scalar []){cc,{-1}})); }
 end_trace("vof_0", "/home/popinet/basilisk-octree/src/vof.h", 203); } return 0; } 
#line 4 "atomisation.c"
#line 1 "tension.h"
#line 1 "/home/popinet/basilisk-octree/src/tension.h"






#line 1 "curvature.h"
#line 1 "/home/popinet/basilisk-octree/src/curvature.h"
#line 12 "/home/popinet/basilisk-octree/src/curvature.h"
static void curvature_restriction (Point point, scalar kappa)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 13 "/home/popinet/basilisk-octree/src/curvature.h"

  double k = 0., s = 0.;
   { foreach_child()
    if (val(kappa,0,0,0) != nodata)
      k += val(kappa,0,0,0), s++; end_foreach_child(); }
  val(kappa,0,0,0) = s ? k/s : nodata;
}







static void curvature_prolongation (Point point, scalar kappa)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 28 "/home/popinet/basilisk-octree/src/curvature.h"

   { foreach_child() {
    double sk = 0., s = 0.;
    for (int i = 0; i <= 1; i++)

      for (int j = 0; j <= 1; j++)


 for (int k = 0; k <= 1; k++)

   if (coarse(kappa,child.x*i,child.y*j,child.z*k) != nodata)
     sk += coarse(kappa,child.x*i,child.y*j,child.z*k), s++;
    val(kappa,0,0,0) = s ? sk/s : nodata;
  } end_foreach_child(); }
}
#line 62 "/home/popinet/basilisk-octree/src/curvature.h"
#line 1 "heights.h"
#line 1 "/home/popinet/basilisk-octree/src/heights.h"
#line 29 "/home/popinet/basilisk-octree/src/heights.h"
static inline double height (double H) {
  return H > 20./2. ? H - 20. : H < -20./2. ? H + 20. : H;
}

static inline int orientation (double H) {
  return fabs(H) > 20./2.;
}
#line 49 "/home/popinet/basilisk-octree/src/heights.h"
static void half_column (Point point, scalar c, vector h, vector cs, int j)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 50 "/home/popinet/basilisk-octree/src/heights.h"







  const int complete = -1;

  {
#line 59
 {







    double S = val(c,0,0,0), H = S, ci, a;







    typedef struct { int s; double h; } HState;
    HState state = {0, 0};
    if (j == 1) {




      if (val(h.x,0,0,0) == 300.)
 state.s = complete, state.h = nodata;




      else {
 int s = (val(h.x,0,0,0) + 20./2.)/100.;
 state.h = val(h.x,0,0,0) - 100.*s;
 state.s = s - 1;
      }





      if (state.s != complete)
 S = state.s, H = state.h;
    }
#line 109 "/home/popinet/basilisk-octree/src/heights.h"
    for (int i = 1; i <= 4; i++) {
      ci = i <= 2 ? val(c,i*j,0,0) : val(cs.x,(i - 2)*j,0,0);
      H += ci;




      if (S > 0. && S < 1.) {
 S = ci;
 if (ci <= 0. || ci >= 1.) {







   H -= i*ci;
   break;
 }
      }
#line 138 "/home/popinet/basilisk-octree/src/heights.h"
      else if (S >= 1. && ci <= 0.) {
 H = (H - 0.5)*j + (j == -1)*20.;
 S = complete;
 break;
      }
      else if (S <= 0. && ci >= 1.) {
 H = (i + 0.5 - H)*j + (j == 1)*20.;
 S = complete;
 break;
      }
#line 156 "/home/popinet/basilisk-octree/src/heights.h"
      else if (S == ci && modf(H, &a))
 break;
    }





    if (j == -1) {







      if (S != complete && ((val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) ||
       (S > 0. && S < 1.)))
 val(h.x,0,0,0) = 300.;
      else if (S == complete)
 val(h.x,0,0,0) = H;
      else





 val(h.x,0,0,0) = H + 100.*(1. + (S >= 1.));
    }
    else {
#line 195 "/home/popinet/basilisk-octree/src/heights.h"
      if (state.s != complete ||
   (S == complete && fabs(height(H)) < fabs(height(state.h))))
 state.s = S, state.h = H;





      if (state.s != complete)
 val(h.x,0,0,0) = nodata;
      else
 val(h.x,0,0,0) = (state.h > 1e10 ? nodata : state.h);
    }
  }
#line 59
 {







    double S = val(c,0,0,0), H = S, ci, a;







    typedef struct { int s; double h; } HState;
    HState state = {0, 0};
    if (j == 1) {




      if (val(h.y,0,0,0) == 300.)
 state.s = complete, state.h = nodata;




      else {
 int s = (val(h.y,0,0,0) + 20./2.)/100.;
 state.h = val(h.y,0,0,0) - 100.*s;
 state.s = s - 1;
      }





      if (state.s != complete)
 S = state.s, H = state.h;
    }
#line 109 "/home/popinet/basilisk-octree/src/heights.h"
    for (int i = 1; i <= 4; i++) {
      ci = i <= 2 ? val(c,0,i*j,0) : val(cs.y,0,(i - 2)*j,0);
      H += ci;




      if (S > 0. && S < 1.) {
 S = ci;
 if (ci <= 0. || ci >= 1.) {







   H -= i*ci;
   break;
 }
      }
#line 138 "/home/popinet/basilisk-octree/src/heights.h"
      else if (S >= 1. && ci <= 0.) {
 H = (H - 0.5)*j + (j == -1)*20.;
 S = complete;
 break;
      }
      else if (S <= 0. && ci >= 1.) {
 H = (i + 0.5 - H)*j + (j == 1)*20.;
 S = complete;
 break;
      }
#line 156 "/home/popinet/basilisk-octree/src/heights.h"
      else if (S == ci && modf(H, &a))
 break;
    }





    if (j == -1) {







      if (S != complete && ((val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) ||
       (S > 0. && S < 1.)))
 val(h.y,0,0,0) = 300.;
      else if (S == complete)
 val(h.y,0,0,0) = H;
      else





 val(h.y,0,0,0) = H + 100.*(1. + (S >= 1.));
    }
    else {
#line 195 "/home/popinet/basilisk-octree/src/heights.h"
      if (state.s != complete ||
   (S == complete && fabs(height(H)) < fabs(height(state.h))))
 state.s = S, state.h = H;





      if (state.s != complete)
 val(h.y,0,0,0) = nodata;
      else
 val(h.y,0,0,0) = (state.h > 1e10 ? nodata : state.h);
    }
  }
#line 59
 {







    double S = val(c,0,0,0), H = S, ci, a;







    typedef struct { int s; double h; } HState;
    HState state = {0, 0};
    if (j == 1) {




      if (val(h.z,0,0,0) == 300.)
 state.s = complete, state.h = nodata;




      else {
 int s = (val(h.z,0,0,0) + 20./2.)/100.;
 state.h = val(h.z,0,0,0) - 100.*s;
 state.s = s - 1;
      }





      if (state.s != complete)
 S = state.s, H = state.h;
    }
#line 109 "/home/popinet/basilisk-octree/src/heights.h"
    for (int i = 1; i <= 4; i++) {
      ci = i <= 2 ? val(c,0,0,i*j) : val(cs.z,0,0,(i - 2)*j);
      H += ci;




      if (S > 0. && S < 1.) {
 S = ci;
 if (ci <= 0. || ci >= 1.) {







   H -= i*ci;
   break;
 }
      }
#line 138 "/home/popinet/basilisk-octree/src/heights.h"
      else if (S >= 1. && ci <= 0.) {
 H = (H - 0.5)*j + (j == -1)*20.;
 S = complete;
 break;
      }
      else if (S <= 0. && ci >= 1.) {
 H = (i + 0.5 - H)*j + (j == 1)*20.;
 S = complete;
 break;
      }
#line 156 "/home/popinet/basilisk-octree/src/heights.h"
      else if (S == ci && modf(H, &a))
 break;
    }





    if (j == -1) {







      if (S != complete && ((val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) ||
       (S > 0. && S < 1.)))
 val(h.z,0,0,0) = 300.;
      else if (S == complete)
 val(h.z,0,0,0) = H;
      else





 val(h.z,0,0,0) = H + 100.*(1. + (S >= 1.));
    }
    else {
#line 195 "/home/popinet/basilisk-octree/src/heights.h"
      if (state.s != complete ||
   (S == complete && fabs(height(H)) < fabs(height(state.h))))
 state.s = S, state.h = H;





      if (state.s != complete)
 val(h.z,0,0,0) = nodata;
      else
 val(h.z,0,0,0) = (state.h > 1e10 ? nodata : state.h);
    }
  }}
}
#line 222 "/home/popinet/basilisk-octree/src/heights.h"
static void column_propagation (vector h)
{
   { foreach(){

#line 224 "/home/popinet/basilisk-octree/src/heights.h"

    for (int i = -2; i <= 2; i++)
      {
#line 226

 if (fabs(height(val(h.x,i,0,0))) <= 3.5 &&
     fabs(height(val(h.x,i,0,0)) + i) < fabs(height(val(h.x,0,0,0))))
   val(h.x,0,0,0) = val(h.x,i,0,0) + i;
#line 226

 if (fabs(height(val(h.y,0,i,0))) <= 3.5 &&
     fabs(height(val(h.y,0,i,0)) + i) < fabs(height(val(h.y,0,0,0))))
   val(h.y,0,0,0) = val(h.y,0,i,0) + i;
#line 226

 if (fabs(height(val(h.z,0,0,i))) <= 3.5 &&
     fabs(height(val(h.z,0,0,i)) + i) < fabs(height(val(h.z,0,0,0))))
   val(h.z,0,0,0) = val(h.z,0,0,i) + i;}; } end_foreach(); }
  boundary ((scalar *)((vector []){{h.x,h.y,h.z},{{-1},{-1},{-1}}}));
}
#line 291 "/home/popinet/basilisk-octree/src/heights.h"

#line 291

static void refine_h_x (Point point, scalar h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 293 "/home/popinet/basilisk-octree/src/heights.h"





  bool complete = true;
   { foreach_child() {
    for (int i = -2; i <= 2; i++)
      if (allocated(i,0,0) && !(!is_leaf(neighbor(i,0,0)) && !neighbor(i,0,0).neighbors && neighbor(i,0,0).pid >= 0) &&
   fabs(height(val(h,i,0,0))) <= 3.5 &&
   fabs(height(val(h,i,0,0)) + i) < fabs(height(val(h,0,0,0))))
 val(h,0,0,0) = val(h,i,0,0) + i;
    if (val(h,0,0,0) == nodata)
      complete = false;
  } end_foreach_child(); }
  if (complete)
    return;
#line 318 "/home/popinet/basilisk-octree/src/heights.h"
  int ori = orientation(val(h,0,0,0));
#line 331 "/home/popinet/basilisk-octree/src/heights.h"
  double H[3][3], H0 = height(val(h,0,0,0));
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      if (val(h,0,i,j) == nodata || orientation(val(h,0,i,j)) != ori)
 return;
      else
 H[i+1][j+1] = height(val(h,0,i,j)) - H0;

  double h0 =
    2.*H0 + (H[2][2] + H[2][0] + H[0][0] + H[0][2] +
      30.*(H[2][1] + H[0][1] + H[1][0] + H[1][2]))/512.
    + 20.*ori;
  double h1 = (H[2][2] + H[2][0] - H[0][0] - H[0][2] +
        30.*(H[2][1] - H[0][1]))/128.;
  double h2 = (H[2][2] - H[2][0] - H[0][0] + H[0][2] +
        30.*(H[1][2] - H[1][0]))/128.;
  double h3 = (H[0][0] + H[2][2] - H[0][2] - H[2][0])/32.;
   { foreach_child()
    if (val(h,0,0,0) == nodata)
      val(h,0,0,0) = h0 + h1*child.y + h2*child.z + h3*child.y*child.z - child.x/2.; end_foreach_child(); }

}
#line 291

static void refine_h_y (Point point, scalar h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 293 "/home/popinet/basilisk-octree/src/heights.h"





  bool complete = true;
   { foreach_child() {
    for (int i = -2; i <= 2; i++)
      if (allocated(0,i,0) && !(!is_leaf(neighbor(0,i,0)) && !neighbor(0,i,0).neighbors && neighbor(0,i,0).pid >= 0) &&
   fabs(height(val(h,0,i,0))) <= 3.5 &&
   fabs(height(val(h,0,i,0)) + i) < fabs(height(val(h,0,0,0))))
 val(h,0,0,0) = val(h,0,i,0) + i;
    if (val(h,0,0,0) == nodata)
      complete = false;
  } end_foreach_child(); }
  if (complete)
    return;
#line 318 "/home/popinet/basilisk-octree/src/heights.h"
  int ori = orientation(val(h,0,0,0));
#line 331 "/home/popinet/basilisk-octree/src/heights.h"
  double H[3][3], H0 = height(val(h,0,0,0));
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      if (val(h,j,0,i) == nodata || orientation(val(h,j,0,i)) != ori)
 return;
      else
 H[i+1][j+1] = height(val(h,j,0,i)) - H0;

  double h0 =
    2.*H0 + (H[2][2] + H[2][0] + H[0][0] + H[0][2] +
      30.*(H[2][1] + H[0][1] + H[1][0] + H[1][2]))/512.
    + 20.*ori;
  double h1 = (H[2][2] + H[2][0] - H[0][0] - H[0][2] +
        30.*(H[2][1] - H[0][1]))/128.;
  double h2 = (H[2][2] - H[2][0] - H[0][0] + H[0][2] +
        30.*(H[1][2] - H[1][0]))/128.;
  double h3 = (H[0][0] + H[2][2] - H[0][2] - H[2][0])/32.;
   { foreach_child()
    if (val(h,0,0,0) == nodata)
      val(h,0,0,0) = h0 + h1*child.z + h2*child.x + h3*child.z*child.x - child.y/2.; end_foreach_child(); }

}
#line 291

static void refine_h_z (Point point, scalar h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 293 "/home/popinet/basilisk-octree/src/heights.h"





  bool complete = true;
   { foreach_child() {
    for (int i = -2; i <= 2; i++)
      if (allocated(0,0,i) && !(!is_leaf(neighbor(0,0,i)) && !neighbor(0,0,i).neighbors && neighbor(0,0,i).pid >= 0) &&
   fabs(height(val(h,0,0,i))) <= 3.5 &&
   fabs(height(val(h,0,0,i)) + i) < fabs(height(val(h,0,0,0))))
 val(h,0,0,0) = val(h,0,0,i) + i;
    if (val(h,0,0,0) == nodata)
      complete = false;
  } end_foreach_child(); }
  if (complete)
    return;
#line 318 "/home/popinet/basilisk-octree/src/heights.h"
  int ori = orientation(val(h,0,0,0));
#line 331 "/home/popinet/basilisk-octree/src/heights.h"
  double H[3][3], H0 = height(val(h,0,0,0));
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      if (val(h,i,j,0) == nodata || orientation(val(h,i,j,0)) != ori)
 return;
      else
 H[i+1][j+1] = height(val(h,i,j,0)) - H0;

  double h0 =
    2.*H0 + (H[2][2] + H[2][0] + H[0][0] + H[0][2] +
      30.*(H[2][1] + H[0][1] + H[1][0] + H[1][2]))/512.
    + 20.*ori;
  double h1 = (H[2][2] + H[2][0] - H[0][0] - H[0][2] +
        30.*(H[2][1] - H[0][1]))/128.;
  double h2 = (H[2][2] - H[2][0] - H[0][0] + H[0][2] +
        30.*(H[1][2] - H[1][0]))/128.;
  double h3 = (H[0][0] + H[2][2] - H[0][2] - H[2][0])/32.;
   { foreach_child()
    if (val(h,0,0,0) == nodata)
      val(h,0,0,0) = h0 + h1*child.x + h2*child.y + h3*child.x*child.y - child.z/2.; end_foreach_child(); }

}







void heights (scalar c, vector h)
{ trace ("heights", "/home/popinet/basilisk-octree/src/heights.h", 361);
  vector cs= new_vector("cs");
  {
#line 363

    for (int i = 0; i < nboundary; i++)
      _attribute[cs.x.i].boundary[i] = _attribute[c.i].boundary[i];
#line 363

    for (int i = 0; i < nboundary; i++)
      _attribute[cs.y.i].boundary[i] = _attribute[c.i].boundary[i];
#line 363

    for (int i = 0; i < nboundary; i++)
      _attribute[cs.z.i].boundary[i] = _attribute[c.i].boundary[i];}





  restriction (((scalar []){c,{-1}}));
  for (int j = -1; j <= 1; j += 2)




    for (int l = 1; l <= depth(); l++) {




       { foreach_level (l){

#line 382 "/home/popinet/basilisk-octree/src/heights.h"

 {
#line 383

   val(cs.x,0,0,0) = val(c,2*j,0,0);
#line 383

   val(cs.y,0,0,0) = val(c,0,2*j,0);
#line 383

   val(cs.z,0,0,0) = val(c,0,0,2*j);}; } end_foreach_level(); }
#line 394 "/home/popinet/basilisk-octree/src/heights.h"
       { foreach_level (l - 1){

#line 394 "/home/popinet/basilisk-octree/src/heights.h"

 {
#line 395
 {
   val(cs.x,0,0,0) = val(c,j,0,0);
   val(cs.x,j,0,0) = val(c,2*j,0,0);
        }
#line 395
 {
   val(cs.y,0,0,0) = val(c,0,j,0);
   val(cs.y,0,j,0) = val(c,0,2*j,0);
        }
#line 395
 {
   val(cs.z,0,0,0) = val(c,0,0,j);
   val(cs.z,0,0,j) = val(c,0,0,2*j);
        }} } end_foreach_level(); }






       { foreach_halo (prolongation, l - 1){

#line 405 "/home/popinet/basilisk-octree/src/heights.h"

 {
#line 406

   _attribute[c.i].prolongation (point, cs.x);
#line 406

   _attribute[c.i].prolongation (point, cs.y);
#line 406

   _attribute[c.i].prolongation (point, cs.z);}; } end_foreach_halo(); }
      { Boundary ** _i = boundaries, * _b; while ((_b = *_i++)) if (_b->halo_prolongation) _b->halo_prolongation (_b, (scalar *)((vector []){{cs.x,cs.y,cs.z},{{-1},{-1},{-1}}}), l, l); };





       { foreach_level (l){

#line 414 "/home/popinet/basilisk-octree/src/heights.h"

        half_column (point, c, h, cs, j); } end_foreach_level(); }
    }






  {
#line 423
 {
    _attribute[h.x.i].prolongation = no_data;
    _attribute[h.x.i].coarsen = no_coarsen;
  }
#line 423
 {
    _attribute[h.y.i].prolongation = no_data;
    _attribute[h.y.i].coarsen = no_coarsen;
  }
#line 423
 {
    _attribute[h.z.i].prolongation = no_data;
    _attribute[h.z.i].coarsen = no_coarsen;
  }}
  boundary ((scalar *)((vector []){{h.x,h.y,h.z},{{-1},{-1},{-1}}}));






  {
#line 434

    _attribute[h.x.i].prolongation = refine_h_x;
#line 434

    _attribute[h.y.i].prolongation = refine_h_y;
#line 434

    _attribute[h.z.i].prolongation = refine_h_z;}




  column_propagation (h);
 delete (((scalar []){cs.x,cs.y,cs.z,{-1}}));  end_trace("heights", "/home/popinet/basilisk-octree/src/heights.h", 441); }
#line 63 "/home/popinet/basilisk-octree/src/curvature.h"
#line 76 "/home/popinet/basilisk-octree/src/curvature.h"

#line 76

static double kappa_z (Point point, vector h) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 77 "/home/popinet/basilisk-octree/src/curvature.h"

  int ori = orientation(val(h.z,0,0,0));
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      if (val(h.z,i,j,0) == nodata || orientation(val(h.z,i,j,0)) != ori)
 return nodata;
  double hx = (val(h.z,1,0,0) - val(h.z,-1,0,0))/2.;
  double hy = (val(h.z,0,1,0) - val(h.z,0,-1,0))/2.;
  double hxx = (val(h.z,1,0,0) + val(h.z,-1,0,0) - 2.*val(h.z,0,0,0))/Delta;
  double hyy = (val(h.z,0,1,0) + val(h.z,0,-1,0) - 2.*val(h.z,0,0,0))/Delta;
  double hxy = (val(h.z,1,1,0) + val(h.z,-1,-1,0) - val(h.z,1,-1,0) - val(h.z,-1,1,0))/(4.*Delta);
  return (hxx*(1. + sq(hy)) + hyy*(1. + sq(hx)) - 2.*hxy*hx*hy)/
    pow(1. + sq(hx) + sq(hy), 3/2.);
}
#line 76

static double kappa_x (Point point, vector h) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 77 "/home/popinet/basilisk-octree/src/curvature.h"

  int ori = orientation(val(h.x,0,0,0));
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      if (val(h.x,0,i,j) == nodata || orientation(val(h.x,0,i,j)) != ori)
 return nodata;
  double hx = (val(h.x,0,1,0) - val(h.x,0,-1,0))/2.;
  double hy = (val(h.x,0,0,1) - val(h.x,0,0,-1))/2.;
  double hxx = (val(h.x,0,1,0) + val(h.x,0,-1,0) - 2.*val(h.x,0,0,0))/Delta;
  double hyy = (val(h.x,0,0,1) + val(h.x,0,0,-1) - 2.*val(h.x,0,0,0))/Delta;
  double hxy = (val(h.x,0,1,1) + val(h.x,0,-1,-1) - val(h.x,0,1,-1) - val(h.x,0,-1,1))/(4.*Delta);
  return (hxx*(1. + sq(hy)) + hyy*(1. + sq(hx)) - 2.*hxy*hx*hy)/
    pow(1. + sq(hx) + sq(hy), 3/2.);
}
#line 76

static double kappa_y (Point point, vector h) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 77 "/home/popinet/basilisk-octree/src/curvature.h"

  int ori = orientation(val(h.y,0,0,0));
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      if (val(h.y,j,0,i) == nodata || orientation(val(h.y,j,0,i)) != ori)
 return nodata;
  double hx = (val(h.y,0,0,1) - val(h.y,0,0,-1))/2.;
  double hy = (val(h.y,1,0,0) - val(h.y,-1,0,0))/2.;
  double hxx = (val(h.y,0,0,1) + val(h.y,0,0,-1) - 2.*val(h.y,0,0,0))/Delta;
  double hyy = (val(h.y,1,0,0) + val(h.y,-1,0,0) - 2.*val(h.y,0,0,0))/Delta;
  double hxy = (val(h.y,1,0,1) + val(h.y,-1,0,-1) - val(h.y,-1,0,1) - val(h.y,1,0,-1))/(4.*Delta);
  return (hxx*(1. + sq(hy)) + hyy*(1. + sq(hx)) - 2.*hxy*hx*hy)/
    pow(1. + sq(hx) + sq(hy), 3/2.);
}
#line 99 "/home/popinet/basilisk-octree/src/curvature.h"
static double height_curvature (Point point, scalar c, vector h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 100 "/home/popinet/basilisk-octree/src/curvature.h"







  typedef struct {
    double n;
    double (* kappa) (Point, vector);
  } NormKappa;
  struct { NormKappa x, y, z; } n;
  {
#line 112

    n.x.n = val(c,1,0,0) - val(c,-1,0,0), n.x.kappa = kappa_x;
#line 112

    n.y.n = val(c,0,1,0) - val(c,0,-1,0), n.y.kappa = kappa_y;
#line 112

    n.z.n = val(c,0,0,1) - val(c,0,0,-1), n.z.kappa = kappa_z;}
  double (* kappaf) (Point, vector) = NULL; NOT_UNUSED (kappaf);




  if (fabs(n.x.n) < fabs(n.y.n))
    swap (NormKappa, n.x, n.y);

  if (fabs(n.x.n) < fabs(n.z.n))
    swap (NormKappa, n.x, n.z);
  if (fabs(n.y.n) < fabs(n.z.n))
    swap (NormKappa, n.y, n.z);





  double kappa = nodata;
  {
#line 132

    if (kappa == nodata) {
      kappa = n.x.kappa (point, h);
      if (kappa != nodata) {
 kappaf = n.x.kappa;
 if (n.x.n < 0.)
   kappa = - kappa;
      }
    }
#line 132

    if (kappa == nodata) {
      kappa = n.y.kappa (point, h);
      if (kappa != nodata) {
 kappaf = n.y.kappa;
 if (n.y.n < 0.)
   kappa = - kappa;
      }
    }
#line 132

    if (kappa == nodata) {
      kappa = n.z.kappa (point, h);
      if (kappa != nodata) {
 kappaf = n.z.kappa;
 if (n.z.n < 0.)
   kappa = - kappa;
      }
    }}

  if (kappa != nodata) {




    if (fabs(kappa) > 1./Delta)
      kappa = sign(kappa)/Delta;
#line 167 "/home/popinet/basilisk-octree/src/curvature.h"
  }

  return kappa;
}
#line 184 "/home/popinet/basilisk-octree/src/curvature.h"
#line 1 "parabola.h"
#line 1 "/home/popinet/basilisk-octree/src/parabola.h"
#line 1 "utils.h"
#line 2 "/home/popinet/basilisk-octree/src/parabola.h"






typedef struct {
  coord o;




  double t[3][3];



  double ** M, rhs[6], a[6];


} ParabolaFit;

static void normalize (coord * n)
{
  double norm = 0.;
  {
#line 26

    norm += sq(n->x);
#line 26

    norm += sq(n->y);
#line 26

    norm += sq(n->z);}
  norm = sqrt(norm);
  {
#line 29

    n->x /= norm;
#line 29

    n->y /= norm;
#line 29

    n->z /= norm;}
}

static void parabola_fit_init (ParabolaFit * p, coord o, coord m)
{
  {
#line 35

    p->o.x = o.x;
#line 35

    p->o.y = o.y;
#line 35

    p->o.z = o.z;}






  double max;
  coord nx = {0., 0., 0.}, ny, nz;
  int d = 0;

  {
#line 47

    nz.x = m.x;
#line 47

    nz.y = m.y;
#line 47

    nz.z = m.z;}
  normalize (&nz);
  max = sq(nz.x);

  if (sq(nz.y) > max) { max = sq(nz.y); d = 1; }
  if (sq(nz.z) > max) d = 2;
  switch (d) {
  case 0: nx.x = - nz.z/nz.x; nx.z = 1.0; break;
  case 1: nx.y = - nz.z/nz.y; nx.z = 1.0; break;
  case 2: nx.z = - nz.x/nz.z; nx.x = 1.0; break;
  }
  normalize (&nx);


  {
#line 62

    ny.x = nz.y*nx.z - nz.z*nx.y;
#line 62

    ny.y = nz.z*nx.x - nz.x*nx.z;
#line 62

    ny.z = nz.x*nx.y - nz.y*nx.x;}


  p->t[0][0] = nx.x; p->t[0][1] = nx.y; p->t[0][2] = nx.z;
  p->t[1][0] = ny.x; p->t[1][1] = ny.y; p->t[1][2] = ny.z;
  p->t[2][0] = nz.x; p->t[2][1] = nz.y; p->t[2][2] = nz.z;



  int n = 6;


  p->M = matrix_new (n, n, sizeof(double));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      p->M[i][j] = 0.;
    p->rhs[i] = 0.;
  }
}

static void parabola_fit_add (ParabolaFit * p, coord m, double w)
{
#line 95 "/home/popinet/basilisk-octree/src/parabola.h"
  double x1 = m.x - p->o.x, y1 = m.y - p->o.y, z1 = m.z - p->o.z;
  double x = p->t[0][0]*x1 + p->t[0][1]*y1 + p->t[0][2]*z1;
  double y = p->t[1][0]*x1 + p->t[1][1]*y1 + p->t[1][2]*z1;
  double z = p->t[2][0]*x1 + p->t[2][1]*y1 + p->t[2][2]*z1;
#line 108 "/home/popinet/basilisk-octree/src/parabola.h"
  double x2 = x*x, x3 = x2*x, x4 = x3*x;
  double y2 = y*y, y3 = y2*y, y4 = y3*y;
  p->M[0][0] += w*x4; p->M[1][1] += w*y4; p->M[2][2] += w*x2*y2;
  p->M[3][3] += w*x2; p->M[4][4] += w*y2; p->M[5][5] += w;
  p->M[0][2] += w*x3*y; p->M[0][3] += w*x3; p->M[0][4] += w*x2*y;
  p->M[1][2] += w*x*y3; p->M[1][3] += w*x*y2; p->M[1][4] += w*y3;
  p->M[2][5] += w*x*y;
  p->M[3][5] += w*x;
  p->M[4][5] += w*y;
  p->rhs[0] += w*x2*z; p->rhs[1] += w*y2*z; p->rhs[2] += w*x*y*z;
  p->rhs[3] += w*x*z; p->rhs[4] += w*y*z; p->rhs[5] += w*z;


}

static double parabola_fit_solve (ParabolaFit * p)
{
#line 149 "/home/popinet/basilisk-octree/src/parabola.h"
  p->M[0][1] = p->M[2][2]; p->M[0][5] = p->M[3][3];
  p->M[1][5] = p->M[4][4];
  p->M[2][3] = p->M[0][4]; p->M[2][4] = p->M[1][3];
  p->M[3][4] = p->M[2][5];
  for (int i = 1; i < 6; i++)
    for (int j = 0; j < i; j++)
      p->M[i][j] = p->M[j][i];
  double pivmin = matrix_inverse (p->M, 6, 1e-10);
  if (pivmin)
    for (int i = 0; i < 6; i++) {
      p->a[i] = 0.;
      for (int j = 0; j < 6; j++)
 p->a[i] += p->M[i][j]*p->rhs[j];
    }
  else
    for (int i = 0; i < 6; i++)
      p->a[i] = 0.;


  matrix_free (p->M);
  return pivmin;
}

static double parabola_fit_curvature (ParabolaFit * p,
          double kappamax, double * kmax)
{
  double kappa;
#line 186 "/home/popinet/basilisk-octree/src/parabola.h"
  double hxx = 2.*p->a[0], hyy = 2.*p->a[1], hxy = p->a[2];
  double hx = p->a[3], hy = p->a[4];

  double dnm = 1. + sq(hx) + sq(hy);
  kappa = - (hxx*(1. + sq(hy)) + hyy*(1. + sq(hx)) - 2.*hxy*hx*hy)
    /sqrt (dnm*dnm*dnm);
  if (kmax) {
    double kg = (hxx*hyy - hxy*hxy)/(dnm*dnm);
    double a = kappa*kappa/4. - kg;
    *kmax = fabs (kappa/2.);
    if (a >= 0.)
      *kmax += sqrt (a);
  }

  if (fabs (kappa) > kappamax) {
    if (kmax)
      *kmax = kappamax;
    return kappa > 0. ? kappamax : - kappamax;
  }
  return kappa;
}
#line 185 "/home/popinet/basilisk-octree/src/curvature.h"






static int independents (coord * p, int n)
{
  if (n < 2)
    return n;
  int ni = 1;
  for (int j = 1; j < n; j++) {
    bool depends = false;
    for (int i = 0; i < j && !depends; i++) {
      double d2 = 0.;
      {
#line 200

 d2 += sq(p[i].x - p[j].x);
#line 200

 d2 += sq(p[i].y - p[j].y);
#line 200

 d2 += sq(p[i].z - p[j].z);}
      depends = (d2 < sq(0.5));
    }
    ni += !depends;
  }
  return ni;
}






static double height_curvature_fit (Point point, scalar c, vector h)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 215 "/home/popinet/basilisk-octree/src/curvature.h"






  coord ip[3 == 2 ? 6 : 27];
  int n = 0;




  {
#line 227
 {





    int n1 = 0, n2 = 0;






    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 if (val(h.z,i,j,0) != nodata) {
   if (orientation(val(h.z,i,j,0))) n1++; else n2++;
 }

    int ori = (n1 > n2);
#line 258 "/home/popinet/basilisk-octree/src/curvature.h"
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 if (val(h.z,i,j,0) != nodata && orientation(val(h.z,i,j,0)) == ori)
   ip[n].x = i, ip[n].y = j, ip[n++].z = height(val(h.z,i,j,0));

  }
#line 227
 {





    int n1 = 0, n2 = 0;






    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 if (val(h.x,0,i,j) != nodata) {
   if (orientation(val(h.x,0,i,j))) n1++; else n2++;
 }

    int ori = (n1 > n2);
#line 258 "/home/popinet/basilisk-octree/src/curvature.h"
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 if (val(h.x,0,i,j) != nodata && orientation(val(h.x,0,i,j)) == ori)
   ip[n].y = i, ip[n].z = j, ip[n++].x = height(val(h.x,0,i,j));

  }
#line 227
 {





    int n1 = 0, n2 = 0;






    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 if (val(h.y,j,0,i) != nodata) {
   if (orientation(val(h.y,j,0,i))) n1++; else n2++;
 }

    int ori = (n1 > n2);
#line 258 "/home/popinet/basilisk-octree/src/curvature.h"
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
 if (val(h.y,j,0,i) != nodata && orientation(val(h.y,j,0,i)) == ori)
   ip[n].z = i, ip[n].x = j, ip[n++].y = height(val(h.y,j,0,i));

  }}





  if (independents (ip, n) < (3 == 2 ? 3 : 9))
    return nodata;





  coord m = mycs (point, c), fc;
  double alpha = plane_alpha (val(c,0,0,0), m);
  double area = plane_area_center (m, alpha, &fc);
  ParabolaFit fit;
  parabola_fit_init (&fit, fc, m);




  parabola_fit_add (&fit, fc, area*100.);






  for (int i = 0; i < n; i++)
    parabola_fit_add (&fit, ip[i], 1.);
  parabola_fit_solve (&fit);
  double kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;



  return kappa;
}






static double centroids_curvature_fit (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 308 "/home/popinet/basilisk-octree/src/curvature.h"






  coord m = mycs (point, c), fc;
  double alpha = plane_alpha (val(c,0,0,0), m);
  plane_area_center (m, alpha, &fc);
  ParabolaFit fit;
  parabola_fit_init (&fit, fc, m);





  coord r = {x,y,z};
   { foreach_neighbor(1)
    if (val(c,0,0,0) > 0. && val(c,0,0,0) < 1.) {
      coord m = mycs (point, c), fc;
      double alpha = plane_alpha (val(c,0,0,0), m);
      double area = plane_area_center (m, alpha, &fc);
      coord rn = {x,y,z};
      {
#line 331

 fc.x += (rn.x - r.x)/Delta;
#line 331

 fc.y += (rn.y - r.y)/Delta;
#line 331

 fc.z += (rn.z - r.z)/Delta;}
      parabola_fit_add (&fit, fc, area);
    } end_foreach_neighbor(); }
  parabola_fit_solve (&fit);
  double kappa = parabola_fit_curvature (&fit, 2., NULL)/Delta;



  return kappa;
}
#line 354 "/home/popinet/basilisk-octree/src/curvature.h"
static inline bool interfacial (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 355 "/home/popinet/basilisk-octree/src/curvature.h"

  if (val(c,0,0,0) >= 1.) {
    for (int i = -1; i <= 1; i += 2)
      {
#line 358

 if (val(c,i,0,0) <= 0.)
   return true;
#line 358

 if (val(c,0,i,0) <= 0.)
   return true;
#line 358

 if (val(c,0,0,i) <= 0.)
   return true;}
  }
  else if (val(c,0,0,0) <= 0.) {
    for (int i = -1; i <= 1; i += 2)
      {
#line 364

 if (val(c,i,0,0) >= 1.)
   return true;
#line 364

 if (val(c,0,i,0) >= 1.)
   return true;
#line 364

 if (val(c,0,0,i) >= 1.)
   return true;}
  }
  else
    return true;
  return false;
}







typedef struct {
  int h;
  int f;
  int a;
  int c;
} cstats;


cstats curvature (scalar c, scalar kappa)
{ trace ("curvature", "/home/popinet/basilisk-octree/src/curvature.h", 388);
  cstats s = {0,0,0,0};
  vector h= new_vector("h");
  heights (c, h);






  _attribute[kappa.i].prolongation = curvature_prolongation;
  _attribute[kappa.i].coarsen = curvature_restriction;






  scalar k= new_scalar("k");
  scalar_clone (k, kappa);

   { foreach(){

#line 409 "/home/popinet/basilisk-octree/src/curvature.h"
 {




    if (!interfacial (point, c))
      val(k,0,0,0) = nodata;





    else if ((val(k,0,0,0) = height_curvature (point, c, h)) != nodata)
      s.h++;
    else if ((val(k,0,0,0) = height_curvature_fit (point, c, h)) != nodata)
      s.f++;
  } } end_foreach(); }
  boundary (((scalar []){k,{-1}}));

   { foreach(){

#line 428 "/home/popinet/basilisk-octree/src/curvature.h"
 {





    if (val(k,0,0,0) < nodata)
      val(kappa,0,0,0) = val(k,0,0,0);
    else if (interfacial (point, c)) {





      double sk = 0., a = 0.;
       { foreach_neighbor(1)
 if (val(k,0,0,0) < nodata)
   sk += val(k,0,0,0), a++; end_foreach_neighbor(); }
      if (a > 0.)
 val(kappa,0,0,0) = sk/a, s.a++;
      else




 val(kappa,0,0,0) = centroids_curvature_fit (point, c), s.c++;
    }
    else
      val(kappa,0,0,0) = nodata;
  } } end_foreach(); }
  boundary (((scalar []){kappa,{-1}}));

  { cstats _ret =  s; delete (((scalar []){k,h.x,h.y,h.z,{-1}}));  end_trace("curvature", "/home/popinet/basilisk-octree/src/curvature.h", 460);  return _ret; }
 delete (((scalar []){k,h.x,h.y,h.z,{-1}}));  end_trace("curvature", "/home/popinet/basilisk-octree/src/curvature.h", 461); }
#line 8 "/home/popinet/basilisk-octree/src/tension.h"











static int defaults_1_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults_1 (const int i, const double t, Event * _ev) { trace ("defaults_1", "/home/popinet/basilisk-octree/src/tension.h", 19);  {







  if (is_constant(a.x)) {
    a = new_face_vector("a");
     { foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 29
{

#line 29 "/home/popinet/basilisk-octree/src/tension.h"

      val(a.x,0,0,0) = 0.; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 29
{

#line 29 "/home/popinet/basilisk-octree/src/tension.h"

      val(a.y,0,0,0) = 0.; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 29
{

#line 29 "/home/popinet/basilisk-octree/src/tension.h"

      val(a.z,0,0,0) = 0.; }  }}  end_foreach_face_generic()
#line 30
 end_foreach_face(); }
    boundary ((scalar *)((vector []){{a.x,a.y,a.z},{{-1},{-1},{-1}}}));
  }





  if (interfaces) for (scalar c = *interfaces, *_i81 = interfaces; ((scalar *)&c)->i >= 0; c = *++_i81)
    if (_attribute[c.i].sigma && !_attribute[c.i].kappa.i) {
      scalar kappa = new_scalar ("kappa");
       { foreach(){

#line 41 "/home/popinet/basilisk-octree/src/tension.h"

 val(kappa,0,0,0) = 0.; } end_foreach(); }
      boundary (((scalar []){kappa,{-1}}));
      _attribute[c.i].kappa = kappa;
    }
 end_trace("defaults_1", "/home/popinet/basilisk-octree/src/tension.h", 46); } return 0; } 
#line 59 "/home/popinet/basilisk-octree/src/tension.h"
static int stability_1_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int stability_1 (const int i, const double t, Event * _ev) { trace ("stability_1", "/home/popinet/basilisk-octree/src/tension.h", 59);  {





  double amin = HUGE, amax = -HUGE, dmin = HUGE;
   { 
#undef _OMPSTART
#define _OMPSTART double _amin = amin; double _amax = amax; double _dmin = dmin; 
#undef _OMPEND
#define _OMPEND OMP(omp critical) if (_amin < amin) amin = _amin; mpi_all_reduce_double (amin, MPI_MIN); OMP(omp critical) if (_amax > amax) amax = _amax; mpi_all_reduce_double (amax, MPI_MAX); OMP(omp critical) if (_dmin < dmin) dmin = _dmin; mpi_all_reduce_double (dmin, MPI_MIN); 
#line 66

if (!is_constant(alpha.x) && !is_constant(fm.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 66
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 66
{

#line 66 "/home/popinet/basilisk-octree/src/tension.h"
 {
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) > _amax) _amax = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) < _amin) _amin = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 66
{

#line 66 "/home/popinet/basilisk-octree/src/tension.h"
 {
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) > _amax) _amax = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) < _amin) _amin = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 66
{

#line 66 "/home/popinet/basilisk-octree/src/tension.h"
 {
    if (val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0) > _amax) _amax = val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0);
    if (val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0) < _amin) _amin = val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  end_foreach_face_generic()
#line 70
 end_foreach_face(); }
if (is_constant(alpha.x) && !is_constant(fm.x)) {
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 66
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 66
{

#line 66 "/home/popinet/basilisk-octree/src/tension.h"
 {
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) > _amax) _amax = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) < _amin) _amin = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 66
{

#line 66 "/home/popinet/basilisk-octree/src/tension.h"
 {
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) > _amax) _amax = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) < _amin) _amin = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 66
{

#line 66 "/home/popinet/basilisk-octree/src/tension.h"
 {
    if (val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0) > _amax) _amax = val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0);
    if (val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0) < _amin) _amin = val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  end_foreach_face_generic()
#line 70
 end_foreach_face(); }
if (!is_constant(alpha.x) && is_constant(fm.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 66
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 66
{

#line 66 "/home/popinet/basilisk-octree/src/tension.h"
 {
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) > _amax) _amax = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) < _amin) _amin = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 66
{

#line 66 "/home/popinet/basilisk-octree/src/tension.h"
 {
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) > _amax) _amax = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) < _amin) _amin = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 66
{

#line 66 "/home/popinet/basilisk-octree/src/tension.h"
 {
    if (val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0) > _amax) _amax = val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0);
    if (val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0) < _amin) _amin = val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  end_foreach_face_generic()
#line 70
 end_foreach_face(); }
if (is_constant(alpha.x) && is_constant(fm.x)) {
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 66
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 66
{

#line 66 "/home/popinet/basilisk-octree/src/tension.h"
 {
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) > _amax) _amax = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0) < _amin) _amin = val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 66
{

#line 66 "/home/popinet/basilisk-octree/src/tension.h"
 {
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) > _amax) _amax = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0) < _amin) _amin = val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 66
{

#line 66 "/home/popinet/basilisk-octree/src/tension.h"
 {
    if (val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0) > _amax) _amax = val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0);
    if (val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0) < _amin) _amin = val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0);
    if (Delta < _dmin) _dmin = Delta;
  } }  }}  end_foreach_face_generic()
#line 70
 end_foreach_face(); }
#undef _OMPSTART
#undef _OMPEND
#define _OMPSTART
#define _OMPEND
#line 71
 }
  double rhom = (1./amin + 1./amax)/2.;





  if (interfaces) for (scalar c = *interfaces, *_i82 = interfaces; ((scalar *)&c)->i >= 0; c = *++_i82)
    if (_attribute[c.i].sigma) {
      double dt = sqrt (rhom*cube(dmin)/(pi*_attribute[c.i].sigma));
      if (dt < dtmax)
 dtmax = dt;
    }
 end_trace("stability_1", "/home/popinet/basilisk-octree/src/tension.h", 83); } return 0; } 
#line 92 "/home/popinet/basilisk-octree/src/tension.h"
static int acceleration_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int acceleration_0 (const int i, const double t, Event * _ev) { trace ("acceleration_0", "/home/popinet/basilisk-octree/src/tension.h", 92); 
{





  scalar * list = NULL;
  if (interfaces) for (scalar c = *interfaces, *_i83 = interfaces; ((scalar *)&c)->i >= 0; c = *++_i83)
    if (_attribute[c.i].sigma) {
      list = list_add (list, c);






       { foreach(){

#line 109 "/home/popinet/basilisk-octree/src/tension.h"

 val(c,0,0,0) = clamp (val(c,0,0,0), 0, 1); } end_foreach(); }
      boundary (((scalar []){c,{-1}}));





      assert (_attribute[c.i].kappa.i);
      curvature (c, _attribute[c.i].kappa);
    }
#line 129 "/home/popinet/basilisk-octree/src/tension.h"
  if (list) for (scalar c = *list, *_i84 = list; ((scalar *)&c)->i >= 0; c = *++_i84)
    _attribute[c.i].prolongation = _attribute[p.i].prolongation;
  boundary (list);
#line 142 "/home/popinet/basilisk-octree/src/tension.h"
  vector st = a;
   { 
if (!is_constant(alpha.x) && !is_constant(fm.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 143
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 143
{

#line 143 "/home/popinet/basilisk-octree/src/tension.h"

    if (list) for (scalar c = *list, *_i85 = list; ((scalar *)&c)->i >= 0; c = *++_i85)
      if (val(c,0,0,0) != val(c,-1,0,0)) {
 scalar kappa = _attribute[c.i].kappa;
#line 156 "/home/popinet/basilisk-octree/src/tension.h"
 double kf =
   (val(kappa,0,0,0) < nodata && val(kappa,-1,0,0) < nodata) ?
      (val(kappa,0,0,0) + val(kappa,-1,0,0))/2. :
   val(kappa,0,0,0) < nodata ? val(kappa,0,0,0) :
   val(kappa,-1,0,0) < nodata ? val(kappa,-1,0,0) :
   0.;

 val(st.x,0,0,0) += val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*_attribute[c.i].sigma*kf*(val(c,0,0,0) - val(c,-1,0,0))/Delta;
      } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 143
{

#line 143 "/home/popinet/basilisk-octree/src/tension.h"

    if (list) for (scalar c = *list, *_i85 = list; ((scalar *)&c)->i >= 0; c = *++_i85)
      if (val(c,0,0,0) != val(c,0,-1,0)) {
 scalar kappa = _attribute[c.i].kappa;
#line 156 "/home/popinet/basilisk-octree/src/tension.h"
 double kf =
   (val(kappa,0,0,0) < nodata && val(kappa,0,-1,0) < nodata) ?
      (val(kappa,0,0,0) + val(kappa,0,-1,0))/2. :
   val(kappa,0,0,0) < nodata ? val(kappa,0,0,0) :
   val(kappa,0,-1,0) < nodata ? val(kappa,0,-1,0) :
   0.;

 val(st.y,0,0,0) += val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*_attribute[c.i].sigma*kf*(val(c,0,0,0) - val(c,0,-1,0))/Delta;
      } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 143
{

#line 143 "/home/popinet/basilisk-octree/src/tension.h"

    if (list) for (scalar c = *list, *_i85 = list; ((scalar *)&c)->i >= 0; c = *++_i85)
      if (val(c,0,0,0) != val(c,0,0,-1)) {
 scalar kappa = _attribute[c.i].kappa;
#line 156 "/home/popinet/basilisk-octree/src/tension.h"
 double kf =
   (val(kappa,0,0,0) < nodata && val(kappa,0,0,-1) < nodata) ?
      (val(kappa,0,0,0) + val(kappa,0,0,-1))/2. :
   val(kappa,0,0,0) < nodata ? val(kappa,0,0,0) :
   val(kappa,0,0,-1) < nodata ? val(kappa,0,0,-1) :
   0.;

 val(st.z,0,0,0) += val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*_attribute[c.i].sigma*kf*(val(c,0,0,0) - val(c,0,0,-1))/Delta;
      } }  }}  end_foreach_face_generic()
#line 164
 end_foreach_face(); }
if (is_constant(alpha.x) && !is_constant(fm.x)) {
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 143
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 143
{

#line 143 "/home/popinet/basilisk-octree/src/tension.h"

    if (list) for (scalar c = *list, *_i85 = list; ((scalar *)&c)->i >= 0; c = *++_i85)
      if (val(c,0,0,0) != val(c,-1,0,0)) {
 scalar kappa = _attribute[c.i].kappa;
#line 156 "/home/popinet/basilisk-octree/src/tension.h"
 double kf =
   (val(kappa,0,0,0) < nodata && val(kappa,-1,0,0) < nodata) ?
      (val(kappa,0,0,0) + val(kappa,-1,0,0))/2. :
   val(kappa,0,0,0) < nodata ? val(kappa,0,0,0) :
   val(kappa,-1,0,0) < nodata ? val(kappa,-1,0,0) :
   0.;

 val(st.x,0,0,0) += val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*_attribute[c.i].sigma*kf*(val(c,0,0,0) - val(c,-1,0,0))/Delta;
      } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 143
{

#line 143 "/home/popinet/basilisk-octree/src/tension.h"

    if (list) for (scalar c = *list, *_i85 = list; ((scalar *)&c)->i >= 0; c = *++_i85)
      if (val(c,0,0,0) != val(c,0,-1,0)) {
 scalar kappa = _attribute[c.i].kappa;
#line 156 "/home/popinet/basilisk-octree/src/tension.h"
 double kf =
   (val(kappa,0,0,0) < nodata && val(kappa,0,-1,0) < nodata) ?
      (val(kappa,0,0,0) + val(kappa,0,-1,0))/2. :
   val(kappa,0,0,0) < nodata ? val(kappa,0,0,0) :
   val(kappa,0,-1,0) < nodata ? val(kappa,0,-1,0) :
   0.;

 val(st.y,0,0,0) += val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*_attribute[c.i].sigma*kf*(val(c,0,0,0) - val(c,0,-1,0))/Delta;
      } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 143
{

#line 143 "/home/popinet/basilisk-octree/src/tension.h"

    if (list) for (scalar c = *list, *_i85 = list; ((scalar *)&c)->i >= 0; c = *++_i85)
      if (val(c,0,0,0) != val(c,0,0,-1)) {
 scalar kappa = _attribute[c.i].kappa;
#line 156 "/home/popinet/basilisk-octree/src/tension.h"
 double kf =
   (val(kappa,0,0,0) < nodata && val(kappa,0,0,-1) < nodata) ?
      (val(kappa,0,0,0) + val(kappa,0,0,-1))/2. :
   val(kappa,0,0,0) < nodata ? val(kappa,0,0,0) :
   val(kappa,0,0,-1) < nodata ? val(kappa,0,0,-1) :
   0.;

 val(st.z,0,0,0) += val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*_attribute[c.i].sigma*kf*(val(c,0,0,0) - val(c,0,0,-1))/Delta;
      } }  }}  end_foreach_face_generic()
#line 164
 end_foreach_face(); }
if (!is_constant(alpha.x) && is_constant(fm.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 143
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 143
{

#line 143 "/home/popinet/basilisk-octree/src/tension.h"

    if (list) for (scalar c = *list, *_i85 = list; ((scalar *)&c)->i >= 0; c = *++_i85)
      if (val(c,0,0,0) != val(c,-1,0,0)) {
 scalar kappa = _attribute[c.i].kappa;
#line 156 "/home/popinet/basilisk-octree/src/tension.h"
 double kf =
   (val(kappa,0,0,0) < nodata && val(kappa,-1,0,0) < nodata) ?
      (val(kappa,0,0,0) + val(kappa,-1,0,0))/2. :
   val(kappa,0,0,0) < nodata ? val(kappa,0,0,0) :
   val(kappa,-1,0,0) < nodata ? val(kappa,-1,0,0) :
   0.;

 val(st.x,0,0,0) += val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*_attribute[c.i].sigma*kf*(val(c,0,0,0) - val(c,-1,0,0))/Delta;
      } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 143
{

#line 143 "/home/popinet/basilisk-octree/src/tension.h"

    if (list) for (scalar c = *list, *_i85 = list; ((scalar *)&c)->i >= 0; c = *++_i85)
      if (val(c,0,0,0) != val(c,0,-1,0)) {
 scalar kappa = _attribute[c.i].kappa;
#line 156 "/home/popinet/basilisk-octree/src/tension.h"
 double kf =
   (val(kappa,0,0,0) < nodata && val(kappa,0,-1,0) < nodata) ?
      (val(kappa,0,0,0) + val(kappa,0,-1,0))/2. :
   val(kappa,0,0,0) < nodata ? val(kappa,0,0,0) :
   val(kappa,0,-1,0) < nodata ? val(kappa,0,-1,0) :
   0.;

 val(st.y,0,0,0) += val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*_attribute[c.i].sigma*kf*(val(c,0,0,0) - val(c,0,-1,0))/Delta;
      } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 143
{

#line 143 "/home/popinet/basilisk-octree/src/tension.h"

    if (list) for (scalar c = *list, *_i85 = list; ((scalar *)&c)->i >= 0; c = *++_i85)
      if (val(c,0,0,0) != val(c,0,0,-1)) {
 scalar kappa = _attribute[c.i].kappa;
#line 156 "/home/popinet/basilisk-octree/src/tension.h"
 double kf =
   (val(kappa,0,0,0) < nodata && val(kappa,0,0,-1) < nodata) ?
      (val(kappa,0,0,0) + val(kappa,0,0,-1))/2. :
   val(kappa,0,0,0) < nodata ? val(kappa,0,0,0) :
   val(kappa,0,0,-1) < nodata ? val(kappa,0,0,-1) :
   0.;

 val(st.z,0,0,0) += val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*_attribute[c.i].sigma*kf*(val(c,0,0,0) - val(c,0,0,-1))/Delta;
      } }  }}  end_foreach_face_generic()
#line 164
 end_foreach_face(); }
if (is_constant(alpha.x) && is_constant(fm.x)) {
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 143
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 143
{

#line 143 "/home/popinet/basilisk-octree/src/tension.h"

    if (list) for (scalar c = *list, *_i85 = list; ((scalar *)&c)->i >= 0; c = *++_i85)
      if (val(c,0,0,0) != val(c,-1,0,0)) {
 scalar kappa = _attribute[c.i].kappa;
#line 156 "/home/popinet/basilisk-octree/src/tension.h"
 double kf =
   (val(kappa,0,0,0) < nodata && val(kappa,-1,0,0) < nodata) ?
      (val(kappa,0,0,0) + val(kappa,-1,0,0))/2. :
   val(kappa,0,0,0) < nodata ? val(kappa,0,0,0) :
   val(kappa,-1,0,0) < nodata ? val(kappa,-1,0,0) :
   0.;

 val(st.x,0,0,0) += val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*_attribute[c.i].sigma*kf*(val(c,0,0,0) - val(c,-1,0,0))/Delta;
      } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 143
{

#line 143 "/home/popinet/basilisk-octree/src/tension.h"

    if (list) for (scalar c = *list, *_i85 = list; ((scalar *)&c)->i >= 0; c = *++_i85)
      if (val(c,0,0,0) != val(c,0,-1,0)) {
 scalar kappa = _attribute[c.i].kappa;
#line 156 "/home/popinet/basilisk-octree/src/tension.h"
 double kf =
   (val(kappa,0,0,0) < nodata && val(kappa,0,-1,0) < nodata) ?
      (val(kappa,0,0,0) + val(kappa,0,-1,0))/2. :
   val(kappa,0,0,0) < nodata ? val(kappa,0,0,0) :
   val(kappa,0,-1,0) < nodata ? val(kappa,0,-1,0) :
   0.;

 val(st.y,0,0,0) += val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*_attribute[c.i].sigma*kf*(val(c,0,0,0) - val(c,0,-1,0))/Delta;
      } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 143
{

#line 143 "/home/popinet/basilisk-octree/src/tension.h"

    if (list) for (scalar c = *list, *_i85 = list; ((scalar *)&c)->i >= 0; c = *++_i85)
      if (val(c,0,0,0) != val(c,0,0,-1)) {
 scalar kappa = _attribute[c.i].kappa;
#line 156 "/home/popinet/basilisk-octree/src/tension.h"
 double kf =
   (val(kappa,0,0,0) < nodata && val(kappa,0,0,-1) < nodata) ?
      (val(kappa,0,0,0) + val(kappa,0,0,-1))/2. :
   val(kappa,0,0,0) < nodata ? val(kappa,0,0,0) :
   val(kappa,0,0,-1) < nodata ? val(kappa,0,0,-1) :
   0.;

 val(st.z,0,0,0) += val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*_attribute[c.i].sigma*kf*(val(c,0,0,0) - val(c,0,0,-1))/Delta;
      } }  }}  end_foreach_face_generic()
#line 164
 end_foreach_face(); } }






  if (list) for (scalar c = *list, *_i86 = list; ((scalar *)&c)->i >= 0; c = *++_i86)
    _attribute[c.i].prolongation = fraction_refine;
  boundary (list);





  pfree (list,__func__,__FILE__,__LINE__);
 end_trace("acceleration_0", "/home/popinet/basilisk-octree/src/tension.h", 180); } return 0; } 
#line 5 "atomisation.c"
#line 21 "atomisation.c"
int maxlevel = 8;
double uemax = 0.1;

scalar f= {11}, * interfaces = ((scalar []){{11},{-1}});
vector alphav= {{12},{13},{14}};
scalar rhov= {15};
vector muv= {{16},{17},{18}};

scalar f0= {19};

static void _set_boundary6 (void) { _attribute[u.x.i].boundary[left] = _boundary6; _attribute[u.x.i].boundary_homogeneous[left] = _boundary6_homogeneous; } 
static void _set_boundary7 (void) { _attribute[u.y.i].boundary[left] = _boundary7; _attribute[u.y.i].boundary_homogeneous[left] = _boundary7_homogeneous; } 
static void _set_boundary8 (void) { _attribute[p.i].boundary[left] = _boundary8; _attribute[p.i].boundary_homogeneous[left] = _boundary8_homogeneous; } 
static void _set_boundary9 (void) { _attribute[f.i].boundary[left] = _boundary9; _attribute[f.i].boundary_homogeneous[left] = _boundary9_homogeneous; } 

static void _set_boundary10 (void) { _attribute[u.x.i].boundary[right] = _boundary10; _attribute[u.x.i].boundary_homogeneous[right] = _boundary10_homogeneous; } 
static void _set_boundary11 (void) { _attribute[p.i].boundary[right] = _boundary11; _attribute[p.i].boundary_homogeneous[right] = _boundary11_homogeneous; } 


static void _set_boundary12 (void) { _attribute[p.i].boundary[top] = _boundary12; _attribute[p.i].boundary_homogeneous[top] = _boundary12_homogeneous; } 
static void _set_boundary13 (void) { _attribute[uf.x.i].boundary[top] = _boundary13; _attribute[uf.x.i].boundary_homogeneous[top] = _boundary13_homogeneous; } 
static void _set_boundary14 (void) { _attribute[p.i].boundary[bottom] = _boundary14; _attribute[p.i].boundary_homogeneous[bottom] = _boundary14_homogeneous; } 
static void _set_boundary15 (void) { _attribute[uf.x.i].boundary[bottom] = _boundary15; _attribute[uf.x.i].boundary_homogeneous[bottom] = _boundary15_homogeneous; } 


static void _set_boundary16 (void) { _attribute[p.i].boundary[front] = _boundary16; _attribute[p.i].boundary_homogeneous[front] = _boundary16_homogeneous; } 
static void _set_boundary17 (void) { _attribute[uf.x.i].boundary[front] = _boundary17; _attribute[uf.x.i].boundary_homogeneous[front] = _boundary17_homogeneous; } 
static void _set_boundary18 (void) { _attribute[p.i].boundary[back] = _boundary18; _attribute[p.i].boundary_homogeneous[back] = _boundary18_homogeneous; } 
static void _set_boundary19 (void) { _attribute[uf.x.i].boundary[back] = _boundary19; _attribute[uf.x.i].boundary_homogeneous[back] = _boundary19_homogeneous; } 




int main (int argc, char * argv[]) { _init_solver();
  if (argc > 1)
    maxlevel = atoi (argv[1]);
  if (argc > 2)
    uemax = atof (argv[2]);

  init_grid (64);

  origin ((struct _origin){0, -1.5, -1.5});
  size (3.);






  alpha = alphav;
  rho = rhov;
  mu = muv;
  _attribute[f.i].sigma = 3e-5;


  mpi.min = 1;
  run();
 free_solver(); }

static int init_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t = 0);   *ip = i; *tp = t;   return ret; } static int init_0 (const int i, const double t, Event * _ev) { trace ("init_0", "atomisation.c", 80);  {

   { foreach(){

#line 82 "atomisation.c"

    val(f,0,0,0) = val(rhov,0,0,0) = val(f0,0,0,0) = 0.; } end_foreach(); }
   { foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 84
{

#line 84 "atomisation.c"

    val(alphav.x,0,0,0) = val(muv.x,0,0,0) = 0.; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 84
{

#line 84 "atomisation.c"

    val(alphav.y,0,0,0) = val(muv.y,0,0,0) = 0.; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 84
{

#line 84 "atomisation.c"

    val(alphav.z,0,0,0) = val(muv.z,0,0,0) = 0.; }  }}  end_foreach_face_generic()
#line 85
 end_foreach_face(); }
  boundary (((scalar []){f,f0,rhov,alphav.x,alphav.y,alphav.z,muv.x,muv.y,muv.z,{-1}}));

  do { int refined; do { refined = 0; ((Quadtree *)grid)->refined.n = 0;  { foreach_leaf(){

#line 88 "atomisation.c"
 if (x < 1.2*0.025 && sq(y) + sq(z) < 2.*sq(1./12.) && level < maxlevel) { refine_cell (point, all, 0, &((Quadtree *)grid)->refined); refined++; continue; } } end_foreach_leaf(); } mpi_all_reduce (refined, MPI_INT, MPI_SUM); if (refined) { mpi_boundary_refine (all, 0); mpi_boundary_update(); boundary (all); while (balance()); } } while (refined); } while(0)
       ;

  scalar phi= new_vertex_scalar("phi");
   { foreach_vertex(){

#line 92 "atomisation.c"

    val(phi,0,0,0) = sq(1./12.) - sq(y) - sq(z); } end_foreach_vertex(); }
  fractions ((struct Fractions){phi, f0});

  _attribute[f0.i].refine = _attribute[f0.i].prolongation = fraction_refine;
  restriction (((scalar []){f0,{-1}}));

   { foreach(){

#line 99 "atomisation.c"
 {
    val(f,0,0,0) = val(f0,0,0,0)*(x < 0.025);
    val(u.x,0,0,0) = val(f,0,0,0);
  } } end_foreach(); }
  boundary (((scalar []){f,u.x,{-1}}));
 delete (((scalar []){phi,{-1}}));  end_trace("init_0", "atomisation.c", 104); } return 0; } 

static int properties_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int properties_0 (const int i, const double t, Event * _ev) { trace ("properties_0", "atomisation.c", 106);  {

  _attribute[f.i].prolongation = refine_bilinear;
  boundary (((scalar []){f,{-1}}));


   { 
if (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 112
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 112
{

#line 112 "atomisation.c"
 {
    double ff = (val(f,0,0,0) + val(f,-1,0,0))/2.;
    val(alphav.x,0,0,0) = val_fm_x(fm.x,0,0,0)/(clamp(ff,0,1)*(1. - 1./27.84) + 1./27.84);
    val(muv.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(2.*1./12./5800*(clamp(ff,0,1)*(1. - 1./27.84) + 1./27.84));
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 112
{

#line 112 "atomisation.c"
 {
    double ff = (val(f,0,0,0) + val(f,0,-1,0))/2.;
    val(alphav.y,0,0,0) = val_fm_y(fm.y,0,0,0)/(clamp(ff,0,1)*(1. - 1./27.84) + 1./27.84);
    val(muv.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(2.*1./12./5800*(clamp(ff,0,1)*(1. - 1./27.84) + 1./27.84));
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 112
{

#line 112 "atomisation.c"
 {
    double ff = (val(f,0,0,0) + val(f,0,0,-1))/2.;
    val(alphav.z,0,0,0) = val_fm_z(fm.z,0,0,0)/(clamp(ff,0,1)*(1. - 1./27.84) + 1./27.84);
    val(muv.z,0,0,0) = val_fm_z(fm.z,0,0,0)*(2.*1./12./5800*(clamp(ff,0,1)*(1. - 1./27.84) + 1./27.84));
  } }  }}  end_foreach_face_generic()
#line 116
 end_foreach_face(); }
if (is_constant(fm.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 112
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 112
{

#line 112 "atomisation.c"
 {
    double ff = (val(f,0,0,0) + val(f,-1,0,0))/2.;
    val(alphav.x,0,0,0) = val_fm_x(fm.x,0,0,0)/(clamp(ff,0,1)*(1. - 1./27.84) + 1./27.84);
    val(muv.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(2.*1./12./5800*(clamp(ff,0,1)*(1. - 1./27.84) + 1./27.84));
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 112
{

#line 112 "atomisation.c"
 {
    double ff = (val(f,0,0,0) + val(f,0,-1,0))/2.;
    val(alphav.y,0,0,0) = val_fm_y(fm.y,0,0,0)/(clamp(ff,0,1)*(1. - 1./27.84) + 1./27.84);
    val(muv.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(2.*1./12./5800*(clamp(ff,0,1)*(1. - 1./27.84) + 1./27.84));
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 112
{

#line 112 "atomisation.c"
 {
    double ff = (val(f,0,0,0) + val(f,0,0,-1))/2.;
    val(alphav.z,0,0,0) = val_fm_z(fm.z,0,0,0)/(clamp(ff,0,1)*(1. - 1./27.84) + 1./27.84);
    val(muv.z,0,0,0) = val_fm_z(fm.z,0,0,0)*(2.*1./12./5800*(clamp(ff,0,1)*(1. - 1./27.84) + 1./27.84));
  } }  }}  end_foreach_face_generic()
#line 116
 end_foreach_face(); } }
   { 
if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 117
foreach(){

#line 117 "atomisation.c"

    val(rhov,0,0,0) = val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0,1)*(1. - 1./27.84) + 1./27.84); } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 117
foreach(){

#line 117 "atomisation.c"

    val(rhov,0,0,0) = val_cm(cm,0,0,0)*(clamp(val(f,0,0,0),0,1)*(1. - 1./27.84) + 1./27.84); } end_foreach(); } }


  _attribute[f.i].prolongation = fraction_refine;
  boundary (((scalar []){f,{-1}}));

 end_trace("properties_0", "atomisation.c", 124); } return 0; } 





static int logfile_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int logfile (const int i, const double t, Event * _ev) { trace ("logfile", "atomisation.c", 130);  {
  if (i == 0)
    fprintf (ferr,
      "t dt mgp.i mgpf.i mgu.i perf.tn perf.t perf.speed mpi.mpe\n");
  fprintf (ferr, "%g %g %d %d %d %ld %g %g %d\n",
    t, dt, mgp.i, mgpf.i, mgu.i,
    perf.tn, perf.t, perf.speed, mpi.npe);
 end_trace("logfile", "atomisation.c", 137); } return 0; } 
#line 160 "atomisation.c"
static int snapshot_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t = 0.1);   *ip = i; *tp = t;   return ret; } static int snapshot_expr1 (int * ip, double * tp, Event * _ev) {   int i = *ip; double t = *tp;   int ret = ( t += 0.1);   *ip = i; *tp = t;   return ret; } static int snapshot_expr2 (int * ip, double * tp, Event * _ev) {   int i = *ip; double t = *tp;   int ret = ( t <= 1.8);   *ip = i; *tp = t;   return ret; } static int snapshot (const int i, const double t, Event * _ev) { trace ("snapshot", "atomisation.c", 160);  {

  char name[80];
  sprintf (name, "snapshot-%g.gfs", t);
  scalar pid= new_scalar("pid"), ff= new_scalar("ff");
   { foreach(){

#line 165 "atomisation.c"
 {
    val(pid,0,0,0) = fmod(pid()*(npe() + 37), npe());
    val(ff,0,0,0) = val(f,0,0,0) < 1e-4 ? 0 : val(f,0,0,0) > 1. - 1e-4 ? 1. : val(f,0,0,0);
  } } end_foreach(); }
  boundary (((scalar []){pid,ff,{-1}}));
  output_gfs ((struct OutputGfs){.file = name, .t = t});

 delete (((scalar []){ff,pid,{-1}}));  end_trace("snapshot", "atomisation.c", 172); } return 0; } 
#line 212 "atomisation.c"
static int adapt_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int adapt (const int i, const double t, Event * _ev) { trace ("adapt", "atomisation.c", 212);  {
#line 225 "atomisation.c"
  boundary ((scalar *)((vector []){{a.x,a.y,a.z},{{-1},{-1},{-1}}}));
  adapt_wavelet ((struct Adapt){((scalar []){f,u.x,u.y,u.z,{-1}}), (double[]){0.01,uemax,uemax,uemax}, maxlevel});
  restriction (((scalar []){f0,{-1}}));
  event ("properties");
#line 246 "atomisation.c"
 end_trace("adapt", "atomisation.c", 246); } return 0; } 
#line 78 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
static double _boundary0 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 77 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 78
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 78
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 78
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 78
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 78
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 78
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 78
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 78
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); } return 0.; } static double _boundary0_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 77 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 78
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 78
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 78
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 78
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 78
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 78
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 78
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 78
return  neumann_homogeneous(); } return 0.; }
#line 79 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
static double _boundary1 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 78 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 79
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 79
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 79
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 79
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); } return 0.; } static double _boundary1_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 78 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 79
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 79
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 79
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 79
return  neumann_homogeneous(); } return 0.; }
#line 85 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
static double _boundary2 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 84 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 85
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 85
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 85
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 85
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 85
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 85
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 85
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 85
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); } return 0.; } static double _boundary2_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 84 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 85
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 85
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 85
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 85
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 85
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 85
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 85
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 85
return  neumann_homogeneous(); } return 0.; }
#line 86 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
static double _boundary3 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 85 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 86
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 86
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 86
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 86
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 86
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 86
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 86
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 86
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); } return 0.; } static double _boundary3_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 85 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 86
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 86
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 86
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 86
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 86
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 86
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 86
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 86
return  neumann_homogeneous(); } return 0.; }
#line 89 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
static double _boundary4 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 88 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); } return 0.; } static double _boundary4_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 88 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann_homogeneous(); } return 0.; }
#line 90 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"
static double _boundary5 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 89 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 90
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 90
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 90
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 90
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 90
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 90
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 90
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 90
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); } return 0.; } static double _boundary5_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 89 "/home/popinet/basilisk-octree/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 90
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 90
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 90
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 90
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 90
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 90
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 90
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 90
return  neumann_homogeneous(); } return 0.; }
#line 31 "atomisation.c"
static double _boundary6 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 30 "atomisation.c"
return  dirichlet(val(f0,0,0,0)*(1. + 0.05*sin (10.*2.*pi*t))); return 0.; } static double _boundary6_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 30 "atomisation.c"
return  dirichlet_homogeneous(); return 0.; }
#line 32 "atomisation.c"
static double _boundary7 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 31 "atomisation.c"
return  dirichlet(0); return 0.; } static double _boundary7_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 31 "atomisation.c"
return  dirichlet_homogeneous(); return 0.; }
#line 33 "atomisation.c"
static double _boundary8 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 32 "atomisation.c"
return  neumann(0); return 0.; } static double _boundary8_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 32 "atomisation.c"
return  neumann_homogeneous(); return 0.; }
#line 34 "atomisation.c"
static double _boundary9 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 33 "atomisation.c"
return  val(f0,0,0,0); return 0.; } static double _boundary9_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 33 "atomisation.c"
return  val(f0,0,0,0); return 0.; }
#line 36 "atomisation.c"
static double _boundary10 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 35 "atomisation.c"
return  neumann(0); return 0.; } static double _boundary10_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 35 "atomisation.c"
return  neumann_homogeneous(); return 0.; }
#line 37 "atomisation.c"
static double _boundary11 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 36 "atomisation.c"
return  dirichlet(0); return 0.; } static double _boundary11_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 36 "atomisation.c"
return  dirichlet_homogeneous(); return 0.; }
#line 40 "atomisation.c"
static double _boundary12 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 39 "atomisation.c"
return  neumann(0); return 0.; } static double _boundary12_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 39 "atomisation.c"
return  neumann_homogeneous(); return 0.; }
#line 41 "atomisation.c"
static double _boundary13 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 40 "atomisation.c"
return  0; return 0.; } static double _boundary13_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 40 "atomisation.c"
return  0; return 0.; }
#line 42 "atomisation.c"
static double _boundary14 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 41 "atomisation.c"
return  neumann(0); return 0.; } static double _boundary14_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 41 "atomisation.c"
return  neumann_homogeneous(); return 0.; }
#line 43 "atomisation.c"
static double _boundary15 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 42 "atomisation.c"
return  0; return 0.; } static double _boundary15_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 42 "atomisation.c"
return  0; return 0.; }
#line 46 "atomisation.c"
static double _boundary16 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 45 "atomisation.c"
return  neumann(0); return 0.; } static double _boundary16_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 45 "atomisation.c"
return  neumann_homogeneous(); return 0.; }
#line 47 "atomisation.c"
static double _boundary17 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 46 "atomisation.c"
return  0; return 0.; } static double _boundary17_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 46 "atomisation.c"
return  0; return 0.; }
#line 48 "atomisation.c"
static double _boundary18 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 47 "atomisation.c"
return  neumann(0); return 0.; } static double _boundary18_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 47 "atomisation.c"
return  neumann_homogeneous(); return 0.; }
#line 49 "atomisation.c"
static double _boundary19 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 48 "atomisation.c"
return  0; return 0.; } static double _boundary19_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 48 "atomisation.c"
return  0; return 0.; }
size_t datasize = 20*sizeof (double);
static int defaults (const int i, const double t, Event * _ev);
static int defaults_expr0 (int * ip, double * tp, Event * _ev);
static int init (const int i, const double t, Event * _ev);
static int init_expr0 (int * ip, double * tp, Event * _ev);
static int set_dtmax (const int i, const double t, Event * _ev);
static int set_dtmax_expr0 (int * ip, double * tp, Event * _ev);
static int stability (const int i, const double t, Event * _ev);
static int stability_expr0 (int * ip, double * tp, Event * _ev);
static int vof (const int i, const double t, Event * _ev);
static int vof_expr0 (int * ip, double * tp, Event * _ev);
static int tracer_advection (const int i, const double t, Event * _ev);
static int tracer_advection_expr0 (int * ip, double * tp, Event * _ev);
static int properties (const int i, const double t, Event * _ev);
static int properties_expr0 (int * ip, double * tp, Event * _ev);
static int advection_term (const int i, const double t, Event * _ev);
static int advection_term_expr0 (int * ip, double * tp, Event * _ev);
static int viscous_term (const int i, const double t, Event * _ev);
static int viscous_term_expr0 (int * ip, double * tp, Event * _ev);
static int acceleration (const int i, const double t, Event * _ev);
static int acceleration_expr0 (int * ip, double * tp, Event * _ev);
static int projection (const int i, const double t, Event * _ev);
static int projection_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_0 (const int i, const double t, Event * _ev);
static int defaults_0_expr0 (int * ip, double * tp, Event * _ev);
static int stability_0 (const int i, const double t, Event * _ev);
static int stability_0_expr0 (int * ip, double * tp, Event * _ev);
static int vof_0 (const int i, const double t, Event * _ev);
static int vof_0_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_1 (const int i, const double t, Event * _ev);
static int defaults_1_expr0 (int * ip, double * tp, Event * _ev);
static int stability_1 (const int i, const double t, Event * _ev);
static int stability_1_expr0 (int * ip, double * tp, Event * _ev);
static int acceleration_0 (const int i, const double t, Event * _ev);
static int acceleration_0_expr0 (int * ip, double * tp, Event * _ev);
static int init_0 (const int i, const double t, Event * _ev);
static int init_0_expr0 (int * ip, double * tp, Event * _ev);
static int properties_0 (const int i, const double t, Event * _ev);
static int properties_0_expr0 (int * ip, double * tp, Event * _ev);
static int logfile (const int i, const double t, Event * _ev);
static int logfile_expr0 (int * ip, double * tp, Event * _ev);
static int snapshot (const int i, const double t, Event * _ev);
static int snapshot_expr0 (int * ip, double * tp, Event * _ev);
static int snapshot_expr1 (int * ip, double * tp, Event * _ev);
static int snapshot_expr2 (int * ip, double * tp, Event * _ev);
static int adapt (const int i, const double t, Event * _ev);
static int adapt_expr0 (int * ip, double * tp, Event * _ev);
static void _set_boundary0 (void);
static void _set_boundary1 (void);
static void _set_boundary2 (void);
static void _set_boundary3 (void);
static void _set_boundary4 (void);
static void _set_boundary5 (void);
static void _set_boundary6 (void);
static void _set_boundary7 (void);
static void _set_boundary8 (void);
static void _set_boundary9 (void);
static void _set_boundary10 (void);
static void _set_boundary11 (void);
static void _set_boundary12 (void);
static void _set_boundary13 (void);
static void _set_boundary14 (void);
static void _set_boundary15 (void);
static void _set_boundary16 (void);
static void _set_boundary17 (void);
static void _set_boundary18 (void);
static void _set_boundary19 (void);
void _init_solver (void) {
  void init_solver();
  init_solver();
  Events = pmalloc (sizeof (Event), __func__, __FILE__, __LINE__);
  Events[0].last = 1;
  event_register ((Event){ 0, 1, defaults, {defaults_expr0}, ((void *)0), ((void *)0),
    "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 97, "defaults"});
  event_register ((Event){ 0, 1, defaults_0, {defaults_0_expr0}, ((void *)0), ((void *)0),
    "/home/popinet/basilisk-octree/src/vof.h", 32, "defaults"});
  event_register ((Event){ 0, 1, defaults_1, {defaults_1_expr0}, ((void *)0), ((void *)0),
    "/home/popinet/basilisk-octree/src/tension.h", 19, "defaults"});
  event_register ((Event){ 0, 1, init, {init_expr0}, ((void *)0), ((void *)0),
    "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 147, "init"});
  event_register ((Event){ 0, 1, init_0, {init_0_expr0}, ((void *)0), ((void *)0),
    "atomisation.c", 80, "init"});
  event_register ((Event){ 0, 1, logfile, {logfile_expr0}, ((void *)0), ((void *)0),
    "atomisation.c", 130, "logfile"});
  event_register ((Event){ 0, 3, snapshot, {snapshot_expr0, snapshot_expr1, snapshot_expr2}, ((void *)0), ((void *)0),
    "atomisation.c", 160, "snapshot"});
  event_register ((Event){ 0, 1, adapt, {adapt_expr0}, ((void *)0), ((void *)0),
    "atomisation.c", 212, "adapt"});
  event_register ((Event){ 0, 1, set_dtmax, {set_dtmax_expr0}, ((void *)0), ((void *)0),
    "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 175, "set_dtmax"});
  event_register ((Event){ 0, 1, stability, {stability_expr0}, ((void *)0), ((void *)0),
    "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 177, "stability"});
  event_register ((Event){ 0, 1, stability_0, {stability_0_expr0}, ((void *)0), ((void *)0),
    "/home/popinet/basilisk-octree/src/vof.h", 44, "stability"});
  event_register ((Event){ 0, 1, stability_1, {stability_1_expr0}, ((void *)0), ((void *)0),
    "/home/popinet/basilisk-octree/src/tension.h", 59, "stability"});
  event_register ((Event){ 0, 1, vof, {vof_expr0}, ((void *)0), ((void *)0),
    "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 187, "vof"});
  event_register ((Event){ 0, 1, vof_0, {vof_0_expr0}, ((void *)0), ((void *)0),
    "/home/popinet/basilisk-octree/src/vof.h", 174, "vof"});
  event_register ((Event){ 0, 1, tracer_advection, {tracer_advection_expr0}, ((void *)0), ((void *)0),
    "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 188, "tracer_advection"});
  event_register ((Event){ 0, 1, properties, {properties_expr0}, ((void *)0), ((void *)0),
    "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 195, "properties"});
  event_register ((Event){ 0, 1, properties_0, {properties_0_expr0}, ((void *)0), ((void *)0),
    "atomisation.c", 106, "properties"});
  event_register ((Event){ 0, 1, advection_term, {advection_term_expr0}, ((void *)0), ((void *)0),
    "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 257, "advection_term"});
  event_register ((Event){ 0, 1, viscous_term, {viscous_term_expr0}, ((void *)0), ((void *)0),
    "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 287, "viscous_term"});
  event_register ((Event){ 0, 1, acceleration, {acceleration_expr0}, ((void *)0), ((void *)0),
    "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 322, "acceleration"});
  event_register ((Event){ 0, 1, acceleration_0, {acceleration_0_expr0}, ((void *)0), ((void *)0),
    "/home/popinet/basilisk-octree/src/tension.h", 92, "acceleration"});
  event_register ((Event){ 0, 1, projection, {projection_expr0}, ((void *)0), ((void *)0),
    "/home/popinet/basilisk-octree/src/navier-stokes/centered.h", 337, "projection"});
  _attribute = pcalloc (datasize/sizeof(double), sizeof (_Attributes), __func__, __FILE__, __LINE__);
  all = pmalloc (sizeof (scalar)*21,__func__, __FILE__, __LINE__);
  for (int i = 0; i < 20; i++)
    all[i].i = i;
  all[20].i = -1;
  set_fpe();
  octree_methods();
  init_scalar ((scalar){19}, "f0");
  init_face_vector ((vector){{16},{17},{18}}, "muv");
  init_scalar ((scalar){15}, "rhov");
  init_face_vector ((vector){{12},{13},{14}}, "alphav");
  init_scalar ((scalar){11}, "f");
  init_face_vector ((vector){{8},{9},{10}}, "uf");
  init_scalar ((scalar){7}, "pf");
  init_vector ((vector){{4},{5},{6}}, "g");
  init_vector ((vector){{1},{2},{3}}, "u");
  init_scalar ((scalar){0}, "p");
  init_const_scalar ((scalar){_NVARMAX+7}, "zeroc",  0.);
  init_const_scalar ((scalar){_NVARMAX+6}, "unity",  1.);
  init_const_vector ((vector){{_NVARMAX+3},{_NVARMAX+4},{_NVARMAX+5}}, "unityf", (double []) {1.,1.,1.});
  init_const_vector ((vector){{_NVARMAX+0},{_NVARMAX+1},{_NVARMAX+2}}, "zerof", (double []) {0.,0.,0.});
  _set_boundary0();
  _set_boundary1();
  _set_boundary2();
  _set_boundary3();
  _set_boundary4();
  _set_boundary5();
  _set_boundary6();
  _set_boundary7();
  _set_boundary8();
  _set_boundary9();
  _set_boundary10();
  _set_boundary11();
  _set_boundary12();
  _set_boundary13();
  _set_boundary14();
  _set_boundary15();
  _set_boundary16();
  _set_boundary17();
  _set_boundary18();
  _set_boundary19();
  init_grid (N);
}
