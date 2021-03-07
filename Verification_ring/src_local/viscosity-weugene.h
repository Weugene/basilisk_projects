#ifndef BASILISK_HEADER_31
#define BASILISK_HEADER_31
#line 1 "./../src_local/viscosity-weugene.h"
#include "penalization.h" // added by Weugene 3/12/20 at 11:55
//#include "poisson.h"
#include "poisson-weugene.h" // added by Weugene 3/12/20 at 11:55

#undef SEPS
#define SEPS 1e-12
bool relative_residual_viscous = false;

#if AXI
#define lambda ((coord){1., 1. + frhs*dt/rho[]*(mu.x[] + mu.x[1] + \
                                mu.y[] + mu.y[0,1])/2./sq(y)})
#else // not AXI
    #if dimension == 1
        #define lambda ((coord){1.})
    #elif dimension == 2
        #define lambda ((coord){1.,1.})
    #elif dimension == 3
        #define lambda ((coord){1.,1.,1.})
    #endif
#endif

#if dimension == 1
    #define scalar_a_by_b(a, b) (a.x[]*b.x[])
#elif dimension == 2
    #define scalar_a_by_b(a, b) (a.x[]*b.x[] + a.y[]*b.y[])
#else // dimension == 3
    #define scalar_a_by_b(a, b) (a.x[]*b.x[] + a.y[]*b.y[] + a.z[]*b.z[])
#endif

#if dimension == 1
    #define m_scalar_a_by_b(a, b) (a.x[]*(b[1] - b[-1]))
#elif dimension == 2
    #define m_scalar_a_by_b(a, b) (a.x[]*(b[1] - b[-1]) + a.y[]*(b[0,1] - b[0,-1]) )
#else // dimension == 3
    #define m_scalar_a_by_b(a, b) (a.x[]*(b[1] - b[-1]) + a.y[]*(b[0,1] - b[0,-1]) + a.z[]*(b[0,0,1] - b[0,0,-1]))
#endif

struct Viscosity {
    vector u;
    face vector mu;
    scalar rho;
    double dt;
    int nrelax;
    scalar * res;
    double maxb;
};

//relaxation function is written for errors e^n = u - u^n, u is an exact solution, u^n is numerical
//here u is e^{n+1}, r is e^n
static void relax_viscosity (scalar*a, scalar*b, int l, void*data)
{
    struct Viscosity * p = (struct Viscosity *) data;
    (const) face vector mu = p->mu;
    (const) scalar rho = p->rho;
    double dt = p->dt;
    vector u = vector(a[0]), r = vector(b[0]);
    coord conv;
    #if JACOBI
        vector w[];
    #else
        vector w = u;
    #endif
//    fprintf(ferr, "relax_viscosity\n");
    foreach_level_or_leaf (l) {
#if BRINKMAN_PENALIZATION
        conv.x = m_scalar_a_by_b(target_U, u.x);
    #if dimension>1
        conv.y = m_scalar_a_by_b(target_U, u.y);
    #endif
    #if dimension>2
        conv.z = m_scalar_a_by_b(target_U, u.z);
    #endif
#endif
        foreach_dimension(){
            w.x[] = (frhs*(dt/rho[])*(2.*mu.x[1]*u.x[1] + 2.*mu.x[]*u.x[-1]
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
                    + r.x[]*sq(Delta))/(sq(Delta)*lambda.x
                    PLUS_DENOMINATOR_BRINKMAN
                    + frhs*(dt/rho[])*(2.*mu.x[1] + 2.*mu.x[]
                #if dimension > 1
                    + mu.y[0,1] + mu.y[]
                #endif
                #if dimension > 2
                    + mu.z[0,0,1] + mu.z[]
                #endif
                    ));
        }
    }

    #if JACOBI
        foreach_level_or_leaf (l)
            foreach_dimension()
                u.x[] = (u.x[] + 2.*w.x[])/3.;
    #endif

    #if TRASH
        vector u1[];
        foreach_level_or_leaf (l)  foreach_dimension() u1.x[] = u.x[];
        trash ({u});
        foreach_level_or_leaf (l) foreach_dimension()  u.x[] = u1.x[];
    #endif
}

//here r is u^n, u is u^{n+1}. u is not an error, it is actual velocity
static double residual_viscosity (scalar * a, scalar * b, scalar * resl, void * data) {
    struct Viscosity * p = (struct Viscosity *) data;
    (const) face vector mu = p->mu;
    (const) scalar rho = p->rho;
    double dt = p->dt;
    vector u = vector(a[0]), r = vector(b[0]), res = vector(resl[0]);
    vector divtauu[];// added: Weugene
    coord conv;// added: Weugene
    double maxres = 0, d = 0, maxb = p->maxb;
    coord LU;
//    fprintf(ferr, "residual_viscosity\n");
    #ifdef CALC_GRAD_U_TAU
        calc_target_U(u, target_U, n_sol);
    #endif
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
        boundary_flux({taux});

        foreach() {
            d = 0.0;
            foreach_dimension() {
                d += taux.x[1] - taux.x[];
            }
            divtauu.x[] = d/Delta;
        }
    }
    boundary((scalar *){divtauu});

    foreach (reduction(max:maxres)) {
#if BRINKMAN_PENALIZATION
        conv.x = m_scalar_a_by_b(target_U, u.x)/(2.0*Delta);
    #if dimension>1
        conv.y = m_scalar_a_by_b(target_U, u.y)/(2.0*Delta);
    #endif
    #if dimension>2
        conv.z = m_scalar_a_by_b(target_U, u.z)/(2.0*Delta);
    #endif
#endif
        /* Lu = b
         * Lu = u^{n+1} - (1-chi)(dt/rho)divtauu + (chi*dt/eta)*u^{n+1} + chi*dt*Ut*\nabla u^{n+1}
         * b = u^n + chi*dt*Ut/eta
         * b is set in viscosity function at the very beginning
         * */
        foreach_dimension(){
            LU.x = lambda.x*u.x[] - frhs*divtauu.x[]*dt/rho[]
                     -(PLUS_VARIABLE_BRINKMAN_RHS)*dt;
            res.x[] = r.x[] - LU.x;
            if (fabs (res.x[]) > maxres) maxres = fabs (res.x[]);
#ifdef DEBUG_BRINKMAN_PENALIZATION
            conv_term.x[] = divtauu.x[]/rho[]; //changed!
            dbp.x[] = PLUS_BRINKMAN_RHS;
            total_rhs.x[] = LU.x;
            residual_of_u.x[] = res.x[];
#endif
        }
    }
    boundary (resl);
    fprintf(ferr, "visc: maxres=%15.12g maxb=%15.12g maxres/maxb=%15.12g\n", maxres, maxb, maxres/maxb);
//    event("vtk_file");
    return maxres/maxb; // Corrected by Weugene: return residual = rhs - du/dt

}

#undef lambda

trace
mgstats viscosity (struct Viscosity p){
    #if AXI
        fprintf(ferr, "Not correct solution for AXI in viscosity-weugene.h");
        exit(9);
    #endif
    #if BRINKMAN_PENALIZATION_ERROR_BC
        fprintf(ferr, "Check BRINKMAN_PENALIZATION_NO_SLIP and BRINKMAN_PENALIZATION_SLIP in viscosity-weugene.h");
        exit(1);
    #endif
    vector u = p.u, r[];
    foreach() foreach_dimension(){ r.x[] = u.x[] + (PLUS_CONSTANT_BRINKMAN_RHS)*dt;}
    face vector mu = p.mu;
    scalar rho = p.rho;
    restriction ({mu,rho});
    if (relative_residual_viscous) {
        double maxb = 0;
        foreach(reduction(max:maxb)){
            if (fabs(r.x[]) > maxb) maxb = fabs(r.x[]);
        }
        if (maxb < SEPS) maxb = 1;
        p.maxb = maxb;
    }else{
        p.maxb = 1;
    }
    fprintf(ferr, "maxb = %g\n", p.maxb);
    return mg_solve ((scalar *){u}, (scalar *){r}, residual_viscosity, relax_viscosity, &p, p.nrelax, p.res);
}

#endif
