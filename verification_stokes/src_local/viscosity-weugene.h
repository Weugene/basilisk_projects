#ifndef BASILISK_HEADER_3
#define BASILISK_HEADER_3
#line 1 "./../src_local/./viscosity-weugene.h"
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
#define fns (1 - fs[])
#else
#define fns 1
#endif


#ifdef BRINKMAN_PENALIZATION
    #if defined(BRINKMAN_PENALIZATION_SLIP) && defined(BRINKMAN_PENALIZATION_NO_SLIP)
        #define  BRINKMAN_PENALIZATION_ERROR_BC
    #endif
    extern scalar fs;
    extern vector Us;
    double eta_s = 1e-3, nu_s = 1e+1;

    #ifdef BRINKMAN_PENALIZATION_NO_SLIP
        #define PLUS_BRINKMAN_RHS         - fs[] * dt * (u.x[] - Us.x[]) / eta_s
        #define PLUS_NUMERATOR_BRINKMAN   + fs[] * dt * sq(Delta) * Us.x[] / eta_s
        #define PLUS_DENOMINATOR_BRINKMAN + fs[] * dt * sq(Delta) / eta_s
    #else
        #define PLUS_BRINKMAN_RHS       - fs[] *  dt * (un * n_sol.x[] / eta_s - nu_s * (u.x[1] - 2 * u.x[] + u.x[-1]) / sq(Delta))
        #if dimension == 1
            #define PLUS_NUMERATOR_BRINKMAN + fs[] * dt * nu_s * (u.x[1] + u.x[-1]) - fs[] * dt * sq(Delta) * (- Usn) * n_sol.x[] / eta_s
        #endif
        #if dimension == 2
            #define PLUS_NUMERATOR_BRINKMAN + fs[] * dt * nu_s * (u.x[1] + u.x[-1]) - fs[] * dt * sq(Delta) * (n_sol.y[] * u.y[] - Usn) * n_sol.x[] / eta_s
        #endif
        #if dimension >= 3
            #define PLUS_NUMERATOR_BRINKMAN + fs[] * dt * nu_s * (u.x[1] + u.x[-1]) - fs[] * dt * sq(Delta) * (n_sol.y[] * u.y[] + n_sol.z[] * u.z[] - Usn) * n_sol.x[] / eta_s
        #endif
        #define PLUS_DENOMINATOR_BRINKMAN   + fs[] * dt * (sq(Delta * n_sol.x[]) / eta_s + 2 * nu_s)
        extern vector n_sol;
    #endif
    #ifdef DEBUG_BRINKMAN_PENALIZATION
        vector dbp[];
    #endif
#endif

static void relax_viscosity (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
  (const) face vector mu = p->mu;
  (const) scalar rho = p->rho;
  double dt = p->dt, Usn = 0;
  vector u = vector(a[0]), r = vector(b[0]);

#if JACOBI
  vector w[];
#else
  vector w = u;
#endif
  
  foreach_level_or_leaf (l) {
    Usn = 0;
    foreach_dimension() Usn += Us.x[]*n_sol.x[];
    foreach_dimension()
      w.x[] = (fns*(dt/rho[])*(2.*mu.x[1]*u.x[1] + 2.*mu.x[]*u.x[-1]
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
                ) PLUS_NUMERATOR_BRINKMAN + r.x[] * sq(Delta))/
               (sq(Delta) * lambda.x + fns * dt / rho[] * (2. * mu.x[1] + 2. * mu.x[]
           #if dimension > 1
				+ mu.y[0,1] + mu.y[]
           #endif
		   #if dimension > 2
				+ mu.z[0,0,1] + mu.z[]
           #endif
                ) PLUS_DENOMINATOR_BRINKMAN);
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
				  void * data) {
    struct Viscosity *p = (struct Viscosity *) data;
    (const) face
    vector mu = p->mu;
    (const) scalar rho = p->rho;
    double dt = p->dt;
    vector u = vector(a[0]), r = vector(b[0]), res = vector(resl[0]);
    double maxres = 0, un = 0, d = 0;
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
#ifdef BRINKMAN_PENALIZATION_SLIP
        un = 0;
        foreach_dimension() un +=(u.x[] - Us.x[]) * n_sol.x[];
#endif
        d = 0;
        foreach_dimension() d += taux.x[1] - taux.x[];
        res.x[] = r.x[]
                  PLUS_BRINKMAN_RHS
                 - lambda.x*u.x[] + fns*dt*d/(Delta*rho[]);
        if (fabs (res.x[]) > maxres)
      maxres = fabs (res.x[]);
      }
    }
    boundary (resl);
#else
    /* "naive" discretisation (only 1st order on trees) */
    foreach(reduction(max:maxres)){
        foreach_dimension()
        {
#ifdef BRINKMAN_PENALIZATION_SLIP
            un = 0;
            foreach_dimension() un += (u.x[] - Us.x[]) * n_sol.x[];
#endif
            res.x[] = r.x[]
            PLUS_BRINKMAN_RHS
            - lambda.x * u.x[]
            + fns * (dt / rho[]) * (2. * mu.x[1, 0] * (u.x[1] - u.x[])
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
            if (fabs(res.x[]) > maxres)
                maxres = fabs(res.x[]);
        }
    }
#endif
#ifdef DEBUG_BRINKMAN_PENALIZATION
    foreach (){
      #ifdef BRINKMAN_PENALIZATION_SLIP
        un = 0;
        foreach_dimension() un +=(u.x[] - Us.x[]) * n_sol.x[];
      #endif
        foreach_dimension() dbp.x[] = PLUS_BRINKMAN_RHS;
    }
#endif
  return maxres;
}

#undef lambda

trace
mgstats viscosity (struct Viscosity p)
{
#if AXI
    fprintf(ferr, "Not correct solution for AXI in viscosity-weugene.h");
    return 9;
#endif
#if BRINKMAN_PENALIZATION_ERROR_BC
    fprintf(ferr, "Check BRINKMAN_PENALIZATION_NO_SLIP and BRINKMAN_PENALIZATION_SLIP in viscosity-weugene.h");
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

#endif
