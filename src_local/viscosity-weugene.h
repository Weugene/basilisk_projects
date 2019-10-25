#include "poisson.h"

struct Viscosity {
  vector u;
  face vector mu;
  scalar rho;
  double dt;
  int nrelax;
  scalar * res;
};

#if AXI
# define lambda ((coord){1., 1. + dt/rho[]*(mu.x[] + mu.x[1] + \
					    mu.y[] + mu.y[0,1])/2./sq(y)})
#else // not AXI
# if dimension == 1
#   define lambda ((coord){1.})
# elif dimension == 2
#   define lambda ((coord){1.,1.})
# elif dimension == 3
#   define lambda ((coord){1.,1.,1.})
#endif
#endif

#ifdef BRINKMAN_PENALIZATION
//extern scalar fs;
//extern vector Us;
double eta_s = 1e-3;
#ifdef DEBUG_BRINKMAN_PENALIZATION
vector dbp[];
#endif
#endif

static void relax_viscosity (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
  (const) face vector mu = p->mu;
  (const) scalar rho = p->rho;
  double dt = p->dt;
  vector u = vector(a[0]), r = vector(b[0]);

#if JACOBI
  vector w[];
#else
  vector w = u;
#endif
  
  foreach_level_or_leaf (l) {
    foreach_dimension()
      w.x[] = (dt/rho[]*(2.*mu.x[1]*u.x[1] + 2.*mu.x[]*u.x[-1]
           #if dimension > 1
			   + mu.y[0,1]*(u.x[0,1] +
					(u.y[1,0] + u.y[1,1])/4. -
					(u.y[-1,0] + u.y[-1,1])/4.)
			   - mu.y[]*(- u.x[0,-1] +
				     (u.y[1,-1] + u.y[1,0])/4. -
				     (u.y[-1,-1] + u.y[-1,0])/4.)
           #endif
	       #if dimension > 2
			   + mu.z[0,0,1]*(u.x[0,0,1] +
					  (u.z[1,0,0] + u.z[1,0,1])/4. -
					  (u.z[-1,0,0] + u.z[-1,0,1])/4.)
			   - mu.z[]*(- u.x[0,0,-1] +
				     (u.z[1,0,-1] + u.z[1,0,0])/4. -
				     (u.z[-1,0,-1] + u.z[-1,0,0])/4.)
           #endif
           #ifdef BRINKMAN_PENALIZATION
			   + fs[]*sq(Delta)*Us.x[]/eta_s
           #endif
               )+ r.x[]*sq(Delta))/
    (sq(Delta)*lambda.x + dt/rho[]*(2.*mu.x[1] + 2.*mu.x[]
           #if dimension > 1
				+ mu.y[0,1] + mu.y[]
           #endif
		   #if dimension > 2
				+ mu.z[0,0,1] + mu.z[]
           #endif
           #ifdef BRINKMAN_PENALIZATION
			    +fs[]*sq(Delta)/eta_s
           #endif
			     ));
  }

#if JACOBI
  foreach_level_or_leaf (l)
    foreach_dimension()
      u.x[] = (u.x[] + 2.*w.x[])/3.;
#endif
  
#if TRASH
  vector u1[];
  foreach_level_or_leaf (l)
    foreach_dimension()
      u1.x[] = u.x[];
  trash ({u});
  foreach_level_or_leaf (l)
    foreach_dimension()
      u.x[] = u1.x[];
#endif
}

static double residual_viscosity (scalar * a, scalar * b, scalar * resl, 
				  void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
  (const) face vector mu = p->mu;
  (const) scalar rho = p->rho;
  double dt = p->dt;
  vector u = vector(a[0]), r = vector(b[0]), res = vector(resl[0]);
  double maxres = 0.;
#if TREE
  /* conservative coarse/fine discretisation (2nd order) */
  foreach_dimension() {
    face vector taux[];
    foreach_face(x)
        taux.x[] = 2.*mu.x[]*(u.x[] - u.x[-1])/Delta;
  #if dimension > 1
    foreach_face(y)
	    taux.y[] = mu.y[]*(u.x[] - u.x[0,-1] +
			   (u.y[1,-1] + u.y[1,0])/4. -
			   (u.y[-1,-1] + u.y[-1,0])/4.)/Delta;
  #endif
  #if dimension > 2
    foreach_face(z)
	    taux.z[] = mu.z[]*(u.x[] - u.x[0,0,-1] +
			   (u.z[1,0,-1] + u.z[1,0,0])/4. -
			   (u.z[-1,0,-1] + u.z[-1,0,0])/4.)/Delta;
  #endif
    boundary_flux ({taux});
    foreach (reduction(max:maxres)) {
      double d = 0.;
      double q = dt/rho[];
    #ifdef BRINKMAN_PENALIZATION
      double qs = fs[]*q/(eta_s);
    #endif
      foreach_dimension()
	    d += taux.x[1] - taux.x[];
      res.x[] = r.x[]
    #ifdef BRINKMAN_PENALIZATION
               + qs*(Us.x[]- u.x[])
    #endif
               - lambda.x*u.x[] + q*d/Delta;
    #ifdef DEBUG_BRINKMAN_PENALIZATION
      dbp.x[] = qs*(Us.x[]- u.x[]);
//      if (fabs(fs[]) > 1e-5) fprintf(stderr, "%g\n", dbp.x[]);
    #endif
      if (fabs (res.x[]) > maxres)
	maxres = fabs (res.x[]);
    }
  }
  boundary (resl);
#else
  /* "naive" discretisation (only 1st order on trees) */
  foreach (reduction(max:maxres)){
        double q = dt/rho[];
    #ifdef BRINKMAN_PENALIZATION
        double qs = fs[]*q/eta_s;
    #endif
        foreach_dimension()
        {
            res.x[] = r.x[]
        #ifdef BRINKMAN_PENALIZATION
                    + qs*(Us.x[]- u.x[])
        #endif
                    - lambda.x * u.x[]
                    + q * (2. * mu.x[1, 0] * (u.x[1] - u.x[])
                         - 2. * mu.x[] * (u.x[] - u.x[-1])
        #if dimension > 1
                    + mu.y[0,1]*(u.x[0,1] - u.x[] +
                               (u.y[1,0] + u.y[1,1])/4. -
                               (u.y[-1,0] + u.y[-1,1])/4.)
                    - mu.y[]*(u.x[] - u.x[0,-1] +
                            (u.y[1,-1] + u.y[1,0])/4. -
                            (u.y[-1,-1] + u.y[-1,0])/4.)
        #endif
        #if dimension > 2
                    + mu.z[0,0,1]*(u.x[0,0,1] - u.x[] +
                             (u.z[1,0,0] + u.z[1,0,1])/4. -
                             (u.z[-1,0,0] + u.z[-1,0,1])/4.)
                    - mu.z[]*(u.x[] - u.x[0,0,-1] +
                            (u.z[1,0,-1] + u.z[1,0,0])/4. -
                            (u.z[-1,0,-1] + u.z[-1,0,0])/4.)
        #endif
                    ) / sq(Delta);
        #ifdef DEBUG_BRINKMAN_PENALIZATION
            dbp.x[] = qs*(Us.x[]- u.x[]);
//            if (fabs(dbp.x[]) > 1e-5) fprintf(stderr, "%g\n", dbp.x[]);
        #endif
            if (fabs(res.x[]) > maxres)
                maxres = fabs(res.x[]);
        }
    }
#endif
  return maxres;
}

#undef lambda

trace
mgstats viscosity (struct Viscosity p)
{
#if AXI
    fprintf(stderr, "Not correct solution for AXI in viscosity-weugene.h");
    return 9;
#endif
  vector u = p.u, r[];
  foreach()
    foreach_dimension()
      r.x[] = u.x[];

  face vector mu = p.mu;
  scalar rho = p.rho;
  restriction ({mu,rho});
  return mg_solve ((scalar *){u}, (scalar *){r},
		   residual_viscosity, relax_viscosity, &p, p.nrelax, p.res);
}

trace
mgstats viscosity_explicit (struct Viscosity p)
{
  vector u = p.u, r[];
  mgstats mg = {0};
  mg.resb = residual_viscosity ((scalar *){u}, (scalar *){u}, (scalar *){r}, &p);
  foreach()
    foreach_dimension()
      u.x[] += r.x[];
  boundary ((scalar *){u});
  return mg;
}
