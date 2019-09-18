/**
# Various utility functions

## Default parameters and variables

The default maximum timestep and CFL number. */

double DT = 1e10, CFL = 0.5;

/**
Performance statistics are stored in this structure. */

struct {
  // total number of (leaf) cells advanced for this process
  long nc;
  // total number of (leaf) cells advanced for all processes
  long tnc;
  // real time elapsed since the start
  double t;
  // average computational speed (leaves/sec)
  double speed;
  // global timer
  timer gt;
} perf;

/**
Performance statistics are gathered by this function, which is
typically called by the run() loop. */

void update_perf() {
  perf.nc += grid->n;
  perf.tnc += grid->tn;
  perf.t = timer_elapsed (perf.gt);
  perf.speed = perf.tnc/perf.t;
}

/**
## Timing

These functions can be used for profiling. */

typedef struct {
  double cpu;   // CPU time (sec)
  double real;  // Wall-clock time (sec)
  double speed; // Speed (points.steps/sec)
  double min;   // Minimum MPI time (sec)
  double avg;   // Average MPI time (sec)
  double max;   // Maximum MPI time (sec)
  size_t tnc;   // Number of grid points
  long   mem;   // Maximum resident memory (kB)
} timing;

/**
Given a timer, iteration count *i*, total number of cells *tnc* and
array of MPI timings *mpi* (with a size equal to the number of
processes), this function returns the statistics above. */

timing timer_timing (timer t, int i, size_t tnc, double * mpi)
{
  timing s;
@if _MPI
  s.avg = mpi_time - t.tm;
@endif
  clock_t end = clock();
  s.cpu = ((double) (end - t.c))/CLOCKS_PER_SEC;
  s.real = timer_elapsed (t);
  if (tnc == 0) {
    double n = 0;
    foreach(reduction(+:n)) n++;
    s.tnc = n;
    tnc = n*i;
  }
  else
    s.tnc = tnc;
@if _GNU_SOURCE
  struct rusage usage;
  getrusage (RUSAGE_SELF, &usage);
  s.mem = usage.ru_maxrss;
@else
  s.mem = 0;
@endif
@if _MPI
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
@else
  s.min = s.max = s.avg = 0.;
@endif
  s.speed = s.real > 0. ? tnc/s.real : -1.;
  return s;
}

/**
This function writes timing statistics on standard output. */

void timer_print (timer t, int i, size_t tnc)
{
  timing s = timer_timing (t, i, tnc, NULL);
  fprintf (fout,
	   "\n# " GRIDNAME 
	   ", %d steps, %g CPU, %.4g real, %.3g points.step/s, %d var\n",
	   i, s.cpu, s.real, s.speed, (int) (datasize/sizeof(double)));
@if _MPI
  fprintf (fout,
	   "# %d procs, MPI: min %.2g (%.2g%%) "
	   "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
	   npe(),
	   s.min, 100.*s.min/s.real,
	   s.avg, 100.*s.avg/s.real,
	   s.max, 100.*s.max/s.real);
@endif
}

/**
## Simple field statistics 

The *normf()* function returns the (volume) average, RMS norm, max
norm and volume for field *f*. */

typedef struct {
  double avg, rms, max, volume;
} norm;

norm normf (scalar f)
{
  double avg = 0., rms = 0., max = 0., volume = 0.;
  foreach(reduction(max:max) reduction(+:avg) 
	  reduction(+:rms) reduction(+:volume)) 
    if (f[] != nodata && dv() > 0.) {
      double v = fabs(f[]);
      if (v > max) max = v;
      volume += dv();
      avg    += dv()*v;
      rms    += dv()*sq(v);
    }
  norm n;
  n.avg = volume ? avg/volume : 0.;
  n.rms = volume ? sqrt(rms/volume) : 0.;
  n.max = max;
  n.volume = volume;
  return n;
}

/**
The *statsf()* function returns the minimum, maximum, volume sum,
standard deviation and volume for field *f*. */

typedef struct {
  double min, max, sum, stddev, volume;
} stats;

stats statsf (scalar f)
{
  double min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
  foreach(reduction(+:sum) reduction(+:sum2) reduction(+:volume)
	  reduction(max:max) reduction(min:min)) 
    if (dv() > 0. && f[] != nodata) {
      volume += dv();
      sum    += dv()*f[];
      sum2   += dv()*sq(f[]);
      if (f[] > max) max = f[];
      if (f[] < min) min = f[];
    }
  stats s;
  s.min = min, s.max = max, s.sum = sum, s.volume = volume;
  if (volume > 0.)
    sum2 -= sum*sum/volume;
  s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
  return s;
}

/**
## Slope limiters 

Given three values, these [slope
limiters](https://en.wikipedia.org/wiki/Flux_limiter#Limiter_functions)
return the corresponding slope-limited gradient. */

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

/**
This is the [generalised minmod
limiter](https://en.wikipedia.org/wiki/Flux_limiter#Generalised_minmod_limiter).
The $\theta$ global variable can be used to tune the limiting
($\theta=1$ gives minmod, the most dissipative limiter and $\theta=2$
gives superbee, the least dissipative). */
	  
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

/**
Given a list of scalar fields *f*, this function fills the gradient
fields *g* with the corresponding gradients. If the *gradient*
attribute of a field is set (typically to one of the limiting
functions above), it is used to compute the gradient, otherwise simple
centered differencing is used. */
	  
void gradients (scalar * f, vector * g)
{
  assert (list_len(f) == vectors_len(g));
  foreach() {
    scalar s; vector v;
    for (s,v in f,g) {
      if (s.gradient)
	foreach_dimension() {
#if EMBED
	  if (!fs.x[] || !fs.x[1])
	    v.x[] = 0.;
	  else
#endif
	    v.x[] = s.gradient (s[-1], s[], s[1])/Delta;
	}
      else // centered
	foreach_dimension() {
#if EMBED
	  if (!fs.x[] || !fs.x[1])
	    v.x[] = 0.;
	  else
#endif
	    v.x[] = (s[1] - s[-1])/(2.*Delta);
	}
    }
  }
  boundary ((scalar *) g);
}

/**
## Other functions

Given a velocity field $\mathbf{u}$, this function fills a scalar
field $\omega$ with the vorticity field
$$
\omega = \partial_x u_y - \partial_y u_x
$$ */

void vorticity (const vector u, scalar omega)
{
  struct { double x, y; } a = {1., -1.};
  foreach() {
    omega[] = 0.;
    foreach_dimension(2)
      omega[] += a.x*center_gradient (u.y);
  }
  boundary ({omega});
}

/**
Given two scalar fields *s* and *sn* this function returns the maximum
of their absolute difference. */

double change (scalar s, scalar sn)
{
  double max = 0.;
  foreach(reduction(max:max)) {
    if (dv() > 0.) {
      double ds = fabs (s[] - sn[]);
      if (ds > max)
	max = ds;
    }
    sn[] = s[];
  }
  return max;
}

#include "output.h"
