#ifndef BASILISK_HEADER_3
#define BASILISK_HEADER_3
#line 1 "./../src_local/./viscosity-weugene.h"
#include "poisson.h"

#if dimension == 1
#define scalar_a_by_b(a, b) (a.x[]*b.x[])
#elif dimension == 2
#define scalar_a_by_b(a, b) (a.x[]*b.x[] + a.y[]*b.y[])
#else // dimension == 3
#define scalar_a_by_b(a, b) (a.x[]*b.x[] + a.y[]*b.y[] + a.z[]*b.z[])
#endif

#if dimension == 1
#define m_scalar_a_by_b(a, b) (a.x[]*(b.x[1] - b.x[-1]))
#elif dimension == 2
#define m_scalar_a_by_b(a, b) (a.x[]*(b.x[1] - b.x[-1]) + a.y[]*(b.y[1] - b.y[-1]))
#else // dimension == 3
#define m_scalar_a_by_b(a, b) (a.x[]*(b.x[1] - b.x[-1]) + a.y[]*(b.y[1] - b.y[-1]) + a.z[]*(b.z[1] - b.z[-1]))
#endif

struct Viscosity {
  vector u;
  face vector mu;
  scalar rho;
  double dt;
  int nrelax;
  scalar*res;
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

    extern scalar fs;
    double eta_s = 1e-3, nu_s = 1e+1;
    (const) vector n_sol = zerof, target_U = zerof;
    (const) scalar a_br = unity;
    (const) scalar b_br = unity;
    vector umUs[];

    #if BRINKMAN_PENALIZATION == 1
        #define PLUS_BRINKMAN_RHS         + fs[]*dt*(nu_s*(u.x[1] - 2*u.x[] +u.x[-1])/sq(Delta) - (a_br[]*u.x[] - target_U.x[])/eta_s)
        #define PLUS_NUMERATOR_BRINKMAN   + fs[]*dt*(nu_s*(u.x[1] + u.x[-1]) + sq(Delta)*target_U.x[]/eta_s)
        #define PLUS_DENOMINATOR_BRINKMAN + fs[]*dt*(a_br[]*sq(Delta) /eta_s + 2*nu_s)
    #elif BRINKMAN_PENALIZATION == 2
        #define CALC_GRAD
        #define PLUS_BRINKMAN_RHS         + fs[]*dt*(nu_s*(u.x[1] - 2*u.x[] +u.x[-1])/sq(Delta) - (b_br[]*scalar_a_by_b(n_sol, grad_u) - target_U.x[])/eta_s)
        #define PLUS_NUMERATOR_BRINKMAN   + fs[]*dt*(nu_s*(u.x[1] + u.x[-1]) - (0.5*Delta*b_br[]*(m_scalar_a_by_b(n_sol, u)) - sq(Delta)*target_U.x[])/eta_s)
        #define PLUS_DENOMINATOR_BRINKMAN + fs[]*dt*(2*nu_s)
    #elif BRINKMAN_PENALIZATION == 3
        #define CALC_GRAD
        #define PLUS_BRINKMAN_RHS         + fs[]*dt*(nu_s*(u.x[1] - 2*u.x[] +u.x[-1])/sq(Delta) - (a_br[]*u.x[] + b_br[]*scalar_a_by_b(n_sol, grad_u) - target_U.x[])/eta_s)
        #define PLUS_NUMERATOR_BRINKMAN   + fs[]*dt*(nu_s*(u.x[1] + u.x[-1]) - (0.5*Delta*b_br[]*(m_scalar_a_by_b(n_sol, u)) - sq(Delta)*target_U.x[])/eta_s)
        #define PLUS_DENOMINATOR_BRINKMAN + fs[]*dt*(a_br[]*sq(Delta)/eta_s + 2*nu_s)
    #elif BRINKMAN_PENALIZATION == 4
        #define CALC_GRAD
        #define PLUS_BRINKMAN_RHS         + fs[]*dt*(nu_s*(u.x[1] - 2*u.x[] +u.x[-1])/sq(Delta) - n_sol.x[]*scalar_a_by_b(n_sol, umUs)/eta_s)
        #define PLUS_NUMERATOR_BRINKMAN   + fs[]*dt*(nu_s*(u.x[1] + u.x[-1]) - sq(Delta)*n_sol.x[]*(scalar_a_by_b(n_sol, umUs) - n_sol.x[]*u.x[])/eta_s)
        #define PLUS_DENOMINATOR_BRINKMAN + fs[]*dt*(sq(n_sol.x[]*Delta)/eta_s + 2*nu_s)
    #else
        #define BRINKMAN_PENALIZATION_ERROR_BC 1
    #endif
    #ifdef DEBUG_BRINKMAN_PENALIZATION
        vector dbp[], total_rhs[];
    #endif
#endif

static void relax_viscosity (scalar*a, scalar*b, int l, void*data)
{
  struct Viscosity*p = (struct Viscosity *) data;
  (const) face vector mu = p->mu;
  (const) scalar rho = p->rho;
  double dt = p->dt;
  vector u = vector(a[0]), r = vector(b[0]);

//    fprintf(ferr, "relax_viscosity...\n");
#if JACOBI
  vector w[];
#else
  vector w = u;
#endif
  
  foreach_level_or_leaf (l) {
    foreach_dimension(){
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
         )
         PLUS_NUMERATOR_BRINKMAN
         + r.x[]*sq(Delta))/
         (sq(Delta)*lambda.x
         PLUS_DENOMINATOR_BRINKMAN
         + fns*dt/rho[]*(2.*mu.x[1] + 2.*mu.x[]
         #if dimension > 1
				 + mu.y[0,1] + mu.y[]
         #endif
		     #if dimension > 2
				 + mu.z[0,0,1] + mu.z[]
         #endif
         ) );
    }
  }
//  foreach_level_or_leaf (l){
//      if (w.x[]>10)fprintf(ferr, "l =%d w= %g u %g %g \n", l, w.x[], u.x[], u.y[]);
//      if (w.x[]>1e15){fprintf(ferr, "l =%d w= %g u %g %g  STOP\n", l, w.x[], u.x[], u.y[]); exit(1);}
//  }
//  fprintf(ferr, "===========================\n");
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

static double residual_viscosity (scalar*a, scalar*b, scalar*resl, 
				  void*data) {
    struct Viscosity *p = (struct Viscosity *) data;
    (const) face
    vector mu = p->mu;
    (const) scalar rho = p->rho;
    double dt = p->dt;
    vector u = vector(a[0]), r = vector(b[0]), res = vector(resl[0]);
    double maxres = 0, un = 0, d = 0;
//    fprintf(ferr, "residual_viscosity...\n");
    /* conservative coarse/fine discretisation (2nd order) */
    foreach_dimension() {
    #ifdef CALC_GRAD
      vector grad_u[];
      gradients({u.x}, {grad_u});
      foreach() umUs.x[] = u.x[] - target_U.x[];
    #endif
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
#ifdef DEBUG_BRINKMAN_PENALIZATION
    foreach_dimension() {
    #ifdef CALC_GRAD
      vector grad_u[];
      gradients({u.x}, {grad_u});
    #endif
      scalar vvv[];

      foreach (){
          dbp.x[] = PLUS_BRINKMAN_RHS;
          total_rhs.x[] = (res.x[] - r.x[] + lambda.x*u.x[])/dt;
//          vvv[]=b_br[]*scalar_a_by_b(n_sol, grad_u);
//          if (dbp.x[]>10) fprintf(ferr, "%g dbp %g rhs %g  grad %g %g \n", vvv[], dbp.x[], total_rhs.x[], grad_u.x[], grad_u.y[]);
//          if (dbp.x[]>10) fprintf(ferr, "b = %g scal_prod %g  \n", b_br[], scalar_a_by_b(n_sol, grad_u));
//          if (dbp.x[]>1e11) {fprintf(ferr, "b = ... big %g scal_prod %g ERROR STOP  \n", b_br[], scalar_a_by_b(n_sol, grad_u));exit(2);}
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





//#else
///* "naive" discretisation (only 1st order on trees) */
//foreach(reduction(max:maxres)){
//foreach_dimension()
//{
//  coord grad_u;
//  foreach_dimension() grad_u.x = (u.x[1] - u.x[-1])/(2.*Delta);
//  res.x[] = r.x[]
//  PLUS_BRINKMAN_RHS
//  - lambda.x*u.x[]
//  + fns*(dt/rho[])*(2.*mu.x[1, 0]*(u.x[1] - u.x[])
//                    - 2.*mu.x[]*(u.x[] - u.x[-1])
//#if dimension > 1
//  + mu.y[0,1]*(u.x[0,1] - u.x[] +
//                       (u.y[1,0] + u.y[1,1])/4. -
//                       (u.y[-1,0] + u.y[-1,1])/4.)
//            - mu.y[]*(u.x[] - u.x[0,-1] +
//                    (u.y[1,-1] + u.y[1,0])/4. -
//                    (u.y[-1,-1] + u.y[-1,0])/4.)
//#endif
//#if dimension > 2
//  + mu.z[0,0,1]*(u.x[0,0,1] - u.x[] +
//                     (u.z[1,0,0] + u.z[1,0,1])/4. -
//                     (u.z[-1,0,0] + u.z[-1,0,1])/4.)
//            - mu.z[]*(u.x[] - u.x[0,0,-1] +
//                    (u.z[1,0,-1] + u.z[1,0,0])/4. -
//                    (u.z[-1,0,-1] + u.z[-1,0,0])/4.)
//#endif
//  )/sq(Delta);
//  if (fabs(res.x[]) > maxres)
//    maxres = fabs(res.x[]);
//}
//}
//#endif
#endif
