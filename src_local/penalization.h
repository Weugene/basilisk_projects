const vector zerocf[] = {0.,0.,0.};
#ifdef BRINKMAN_PENALIZATION
    #define frhs (1)
    #define fbp (fs[])
    extern scalar fs;
    extern face vector fs_face;
    double eta_s = 1e-15, nu_s = 1, lambda_slip = 0;
    double m_bp = 0;
    (const) scalar a_br = unity, b_br = unity; // useful for Robin BC
    (const) vector U_solid = zerocf;
    #if BRINKMAN_PENALIZATION == 1 //Dirichlet BC
        (const) vector target_U = zerocf, n_sol = zerocf;
//        (const) face vector target_Uf = zerof;
        #define PLUS_CONSTANT_BRINKMAN_RHS + (fbp*(target_U.x[])/eta_s)
        #define PLUS_VARIABLE_BRINKMAN_RHS - (fbp*(u.x[]/eta_s))
        #define PLUS_NUMERATOR_BRINKMAN    + 0
        #define PLUS_DENOMINATOR_BRINKMAN  + fbp*dt*(sq(Delta)/eta_s)
    #elif BRINKMAN_PENALIZATION == 2 //Neumann BC
        vector target_U[], n_sol[];
        #define CALC_GRAD
        #define PLUS_CONSTANT_BRINKMAN_RHS  + (fbp*(target_U.x[])/eta_s)
        #define PLUS_VARIABLE_BRINKMAN_RHS  + fbp*dt*( - (scalar_a_by_b(n_sol, grad_u) )/eta_s)
        #define PLUS_NUMERATOR_BRINKMAN   + fbp*dt*( - (0.5*Delta*(m_scalar_a_by_b(n_sol, u)) - sq(Delta)*target_U.x[])/eta_s)
        #define PLUS_DENOMINATOR_BRINKMAN + 0
    #elif BRINKMAN_PENALIZATION == 3 //Robin BC
        vector target_U[], n_sol[];
        #define CALC_GRAD
        #define PLUS_CONSTANT_BRINKMAN_RHS  + (fbp*(target_U.x[])/eta_s)
        #define PLUS_VARIABLE_BRINKMAN_RHS + fbp*dt*( - (a_br[]*u.x[] + b_br[]*scalar_a_by_b(n_sol, grad_u) )/eta_s)
        #define PLUS_NUMERATOR_BRINKMAN   + fbp*dt*( - (0.5*Delta*b_br[]*(m_scalar_a_by_b(n_sol, u)) - sq(Delta)*target_U.x[])/eta_s)
        #define PLUS_DENOMINATOR_BRINKMAN + fbp*dt*(a_br[]*sq(Delta)/eta_s)
    #elif BRINKMAN_PENALIZATION == 4 //No penetration & slip BC
        vector target_U[], n_sol[]; //here target_U will be recalcilated each time step
        #define CALC_GRAD_U_TAU
        #define PLUS_CONSTANT_BRINKMAN_RHS + (fbp*(target_U.x[])/eta_s)
        #define PLUS_VARIABLE_BRINKMAN_RHS - (fbp*(u.x[]/eta_s ))
        #define PLUS_NUMERATOR_BRINKMAN    + 0
        #define PLUS_DENOMINATOR_BRINKMAN  + fbp*dt*(sq(Delta)/eta_s)
//    #elif BRINKMAN_PENALIZATION == 5//Ideal slip, no friction
//        #define CALC_GRAD
//        #define PLUS_BRINKMAN_RHS         + fbp*dt*(nu_s*(u.x[1] - 2*u.x[] +u.x[-1])/sq(Delta) - n_sol.x[]*scalar_a_by_b(n_sol, umUs)/eta_s)
//        #define PLUS_NUMERATOR_BRINKMAN   + fbp*dt*(nu_s*(u.x[1] + u.x[-1]) - sq(Delta)*n_sol.x[]*(scalar_a_by_b(n_sol, umUs) - n_sol.x[]*u.x[])/eta_s)
//        #define PLUS_DENOMINATOR_BRINKMAN + fbp*dt*(sq(n_sol.x[]*Delta)/eta_s + 2*nu_s)
    #else
        #define BRINKMAN_PENALIZATION_ERROR_BC 1
    #endif
    #undef SEPS
    #define SEPS 1e-15
    #ifdef DEBUG_BRINKMAN_PENALIZATION
        vector dbp[], total_rhs[], residual_of_u[], divtauu[];
        vector utau[], grad_utau_n[];
        #define gradun grad_utau_n.x[]
    #else
        double gradun;
    #endif
#else
    #define frhs 1
    #define fbp 0
    double gradun;
    #define PLUS_CONSTANT_BRINKMAN_RHS  0
    #define PLUS_VARIABLE_BRINKMAN_RHS 0
    #define PLUS_NUMERATOR_BRINKMAN 0
    #define PLUS_DENOMINATOR_BRINKMAN 0
#endif

#define PLUS_BRINKMAN_RHS  (PLUS_VARIABLE_BRINKMAN_RHS + PLUS_CONSTANT_BRINKMAN_RHS) //- fbp*( (u.x[] - target_U.x[])/eta_s )

struct Brinkman {
    vector u;
    face vector uf;
    scalar rho;
    double dt;
};

#if BRINKMAN_PENALIZATION == 4
void calc_target_U(const vector u, vector target_U, const vector normal){
    if (!is_constant(U_solid.x)) foreach() foreach_dimension() target_U.x[] = U_solid.x[];
    if (fabs(lambda_slip) > 0.) {
        double ubyn;
        #ifndef DEBUG_BRINKMAN_PENALIZATION
        vector utau[]; // otherwise utau will be defined globally
        #endif
        if (!is_constant(U_solid.x)) foreach() foreach_dimension() u.x[] -= U_solid.x[];
        //    if (!is_constant(U_solid.x)) fprintf(ferr, "U_solid.x");
        foreach() {
            ubyn = 0;
            foreach_dimension() ubyn += u.x[]*normal.x[];
            foreach_dimension() utau.x[] = u.x[] - ubyn*normal.x[];
        }
        if (!is_constant(U_solid.x)) foreach() foreach_dimension() u.x[] += U_solid.x[];

        foreach() {
            if (0 < fs[] && fs[] < 1) { //See here!
				coord gradun;
                foreach_dimension() {
					gradun.x = (n_sol.x[] > 0) ? (utau.x[2] - utau.x[1])/Delta : (utau.x[-1] - utau.x[-2])/Delta;
    //                    gradun = ((1 - fs[-1])*utau.x[-1] - (1 - fs[-2])*utau.x[-2] + (1 - fs[2])*utau.x[2] - (1 - fs[1])*utau.x[1])*n_sol.x[] / Delta
    //                            #if dimension > 1
    //                            + ((1 - fs[0,-1])*utau.x[0,-1] - (1 - fs[0,-2])*utau.x[0,-2] + (1 - fs[0,2])*utau.x[0,2] - (1 - fs[0,1])*utau.x[0,1])*n_sol.y[]/Delta
    //                            #endif
    //                            #if dimension > 2
    //                            + ((1 - fs[0, 0, -1])*utau.x[0, 0, -1] - (1 - fs[0, 0, -2])*utau.x[0, 0, -2] +  (1 - fs[0, 0, 2])*utau.x[0, 0, 2] - (1 - fs[0, 0, 1])*utau.x[0, 0, 1])*n_sol.z[] / Delta
    //                            #endif
    //                            ;
    //                    gradun = (utau.x[-1] - utau.x[-2] + utau.x[2] -  utau.x[1])*n_sol.x[] / Delta
    //                            #if dimension > 1
    //                            + (utau.x[0,-1] - utau.x[0,-2] + utau.x[0,2] - utau.x[0,1])*n_sol.y[]/Delta
    //                            #endif
    //                            #if dimension > 2
    //                            + (utau.x[0, 0, -1] - utau.x[0, 0, -2] +  utau.x[0, 0, 2] - utau.x[0, 0, 1])*n_sol.z[] / Delta
    //                            #endif
    //                            ;
                        gradun = ((utau.x[-1] - utau.x[-2])*0.5*(fs[-1] + fs[-2]) + (utau.x[2] -  utau.x[1])*0.5*(fs[1] + fs[2]))*n_sol.x[] / Delta
                                #if dimension > 1
                                + ((utau.x[0,-1] - utau.x[0,-2])*0.5*(fs[0,-1] + fs[0,-2]) + (utau.x[0,2] - utau.x[0,1])*0.5*(fs[0,1] + fs[0,2]))*n_sol.y[]/Delta
                                #endif
                                #if dimension > 2
                                + ((utau.x[0, 0, -1] - utau.x[0, 0, -2])*0.5*(fs[0, 0, -1] + fs[0, 0, -2]) + (utau.x[0, 0, 2] - utau.x[0, 0, 1])*0.5*(fs[0, 0, 2] + fs[0, 0, 1]))*n_sol.z[] / Delta
                                #endif
                                ;
                        target_U.x[] += lambda_slip*gradun;
                }
                //            if (target_U.x[]) fprintf(ferr, "Ut = %g, U_s = %g, l=%g gr=%g\n", target_U.x[], U_solid.x[], lambda_slip, gradun);
            }else{
                foreach_dimension() target_U.x[] = u.x[];
            }
        }
    }
}
#endif
void brinkman_correction_u (vector u, double dt){
#if BRINKMAN_PENALIZATION == 4
    if (!is_constant(target_U.x)) calc_target_U(u, target_U, n_sol);
#endif
    foreach() {
        foreach_dimension(){
            u.x[] = (u.x[] + (fbp*dt/eta_s)*target_U.x[])/(1. + fbp*dt/eta_s);
        }
    }
    boundary ((scalar *){u});
}

void brinkman_correction_uf (face vector uf){
    foreach_face() {
        uf.x[] = (1 - fs_face.x[])*uf.x[] + fs_face.x[]*face_value(target_U.x, 0);
    }
    boundary((scalar *){uf});
}

void brinkman_correction (struct Brinkman p){
    vector u = p.u; double dt = p.dt;
    brinkman_correction_u (u, dt);
//    brinkman_correction_uf (uf, dt);
}
