//The module is intended to solve the heat transfer equation, the polymerization effect and rheology changing
//it is coupled with navier-stokes/centered.h module
//
#define HEAT_TRANSFER
#define REACTION_MODEL_NON_AUTOCATALYTIC 1
#define REACTION_MODEL_N_ORDER_AUTOCATALYTIC 2
#define REACTION_MODEL_PROUT_TOMPKINS_AUTOCATALYTIC 3
#define NO_REACTION_MODEL 4

#include "three-phase-rheology.h"
#include "diffusion-weugene.h"
#include "dissipation.h"
//const scalar const_temp_solid[] = 1.1;
(const) scalar T_target = unity;
double eta_T = 0;
double m_bp_T = 0;
double chi_conductivity = 0;
double Htr = 1;
double K_cat = 1;
double Arrhenius_const = 10;//1/s
double Ea_by_R = 5;// Kelvin
double n_degree = 1.667;
double m_degree = 0.333;
bool viscDissipation = false;
#undef SEPS
#define SEPS 1e-10

/**
 * In the DSC measurements, the degree of cure ($\alpha$) ranges from 0 (completely uncured) to 1 (fully cured) and is defined as follows:
 * $\alpha(t) = \frac{H(t)}{H_{tr}},$
 * where $H(t)$ is the heat of the reaction up to time t and $H_{tr}$ refers to the total heat of the reaction.
 * Common  Differential Scanning Calorimetry (DSC) measurements involve measurements of the heat flow as a function of time, and
 * the enthalpy can be obtained by integrating the area of exothermic peak,
 * which is the integration of the heat flow.
 * Once the cure behavior of the resin has been examined with the DSC, it is
 * important to find a kinetics model that can describe the cure behavior.
 * A variety of kinetics models have been developed to relate the rate of cure and degree of cure.
 * Phenomenological reaction models are the most common models to describe thermoset cure reactions.
 * The kinetic parameters of the cure reaction can be obtained by fitting the data obtained
 * from the DSC measurements to the phenomenological reaction models.
 * Kinetic models:
 * It is assumed that the rate of reaction can be defined by two separable functions
 * $$\frac{d\alpha}{d\alpha} = K(T) FR(\alpha)$$
 * where $\frac{d\alpha}{d\alpha}$ is the rate of reaction, $K(T)$ is the temperature dependent rate constant,
 * and $FR(\alpha)$ corresponds to the reaction model.
 * The temperature dependence of the reaction rate is generally defined through an Arrhenius expression:
 * $$K(T) = A \exp(\frac{-E_a}{RT}),$$
 * where $E_a$ is the activation energy, $A$ is the pre-exponential factor,
 * $R$ refers to the universal gas constant, and $T$ corresponds to the absolute temperature.

 */

#define KT(T) ( Arrhenius_const*exp(-Ea_by_R/T) )
#define dKT_dT(T) ( KT(T)*Ea_by_R/sq(T) )


#if REACTION_MODEL == REACTION_MODEL_N_ORDER_AUTOCATALYTIC
    #define FR(alpha_doc) ( pow(1 - alpha_doc, n_degree)*(1 + K_cat*alpha_doc) )
    #define dFR_dalpha(alpha_doc) ( pow(1 - alpha_doc, n_degree)*( -n_degree/(1 - alpha_doc) + K_cat) )
    #define GENERAL_METHOD 1
#elif REACTION_MODEL == REACTION_MODEL_PROUT_TOMPKINS_AUTOCATALYTIC
    #define FR(alpha_doc) ( pow(1 - alpha_doc, n_degree)*pow(alpha_doc, m_degree) )
    #define dFR_dalpha(alpha_doc) ( FR(alpha_doc) * ( -n_degree/(1 - alpha_doc) + m_degree/alpha_doc) )
    #define GENERAL_METHOD 1
#else //REACTION_MODEL_NON_AUTOCATALYTIC
    #define FR(alpha_doc) ( pow(1 - alpha_doc, n_degree) )
    #define dFR_dalpha(alpha_doc) ( -n_degree * pow(1 - alpha_doc, n_degree - 1) )
    #define GENERAL_METHOD 0
#endif

/**
 * \rho C_p T_t = \nabla\kappa\nabla T^{n+1} + \rho_1 Q A (1-\alpha^n)^{n_degree}\exp(-E_a/(RT^n))(1 - E_a/(R T^n) + E_a T^{n+1}/(R (T^n)^2)) - \rho C_p \chi\frac{T^{n+1}- T_0}{\eta_T}
 * \thetav = \rho C_p
 * D = \kappav
 * beta = \rho_1 Q A (1-\alpha^n)^{n_degree} \exp(-E_a/(RT^n)) \frac{E_a}{R (T^n)^2} - \frac{\rho C_p \chi}{\eta_T}
 * r = \rho_1 Q A (1-\alpha^n)^{n_degree} \exp(-\frac{E_a}{RT^n})(1 - \frac{E_a}{RT^n}) + \frac{\rho C_p \chi T_0}{\eta_T}
 */


#if REACTION_MODEL != NO_REACTION_MODEL //  POLYMERIZATION_REACTION
event vof (i++){
    if (!stokes_heat) {
        advection ((scalar *){alpha_doc}, uf, dt);
    }
}
#endif

event chem_advection_term (i++){
    if (!stokes_heat) {
        advection((scalar *) {T}, uf, dt);
    }
}

mgstats mgT;
event end_timestep (i++){
    scalar r[], thetav[], beta[], R_source[];
    foreach()
        thetav[] =  var_hom(f[], fs[], rho1*Cp1, rho2*Cp2, rho3*Cp3);
    // advection term of kinetic equation is solved implicitly
    // due to a linearized source term
#if REACTION_MODEL != NO_REACTION_MODEL //  POLYMERIZATION_REACTION
    foreach() {
        double alpha_doc_old = alpha_doc[];
#if GENERAL_METHOD == 0
        alpha_doc[] = 1.0 - pow(
                pow(fabs(1 - alpha_doc[]), 1 - n_degree) + (n_degree - 1.0) * dt * KT(T[]),
                1.0/(1.0 - n_degree));//direct integration from t to t + dt at fixed T
#else
// Crank--Nicolson
//        alpha_doc[] = (alpha_doc[] + 0.5 * dt * KT(T[]) * (2.0 * FR(alpha_doc[]) - dFR_dalpha(alpha_doc[]) * alpha_doc[])) /
//                          (1 - 0.5 * dt * KT(T[]) * dFR_dalpha(alpha_doc[]));
// Backward Euler
        alpha_doc[] = (alpha_doc[] + dt * KT(T[]) * ( FR(alpha_doc[]) - dFR_dalpha(alpha_doc[]) * alpha_doc[])) /
                          (1 - dt * KT(T[]) * dFR_dalpha(alpha_doc[]));
#endif
        //source and conductivity terms. f[] * (1 - fs[]) multiplications means gas and solids can't produce heat
        R_source[] = (alpha_doc[] - alpha_doc_old) * f[] * (1 - fs[]) / dt;
        alpha_doc[] = clamp(alpha_doc[], 0.0, 1.0);
    }
#endif
    // considering the viscous dissipation
    if ((i>1000) && viscDissipation){
        dissipation (r, u, mu = mu);
    }else{
        foreach () {
            r[] = 0.;
        }
    }
    //solids play role in the conduction process
    if (fabs(Htr) > SEPS){
        foreach() {
            beta[] = 0;
            #if REACTION_MODEL != NO_REACTION_MODEL
                r[] += rho1*Htr*R_source[];
            #endif
            //Penalization terms are:
            #if T_DIRICHLET_BC == 1
                r[] += fs[] * thetav[] * T_target[] / eta_T;
                beta[] += -fs[] * thetav[] / eta_T;
            #endif
        }
    }else {// Htr = 0
            foreach() {
                #if T_DIRICHLET_BC == 1
                    r[] += fs[] * thetav[] * T_target[] / eta_T;
                    beta[] = -fs[] * thetav[] / eta_T;
                #else
                    beta[] = 0;
                #endif
            }
    }

    if (constant(kappa.x) != 0.) {
        mgT = diffusion(T, dt, D = kappav, r = r, beta = beta, theta = thetav);
#ifdef DEBUG_HEAT
        fprintf (stderr, "mgT: i=%d t=%g dt=%g num of iterations T=%d\n", i, t, dt, mgT.i); //number of iterations
#endif
    }else{
        foreach(){
            T[] = (thetav[] * T[] + dt * r[])/(thetav[] - dt * beta[]);
        }
    }
}

// In first 10 steps, eta_T and m_bp_T will be adjusted
event properties (i < 10){
    chi_conductivity = kappa1 / (rho1 * Cp1);

    double mindelta=1e+10;
    foreach( reduction(min:mindelta) ){
        if (Delta < mindelta) mindelta = Delta;
    }
	fprintf(ferr, "chi_conductivity=%g\n", chi_conductivity);
    if (chi_conductivity) {
        if (fabs(eta_T) > SEPS && fabs(m_bp_T) > 0) { // m_bp_T has higher priority
            eta_T = sq(m_bp_T * mindelta) / chi_conductivity;
        } else if (fabs(eta_T) < SEPS && fabs(m_bp_T) == 0) { // nothing is set
            m_bp_T = 1;
            eta_T = sq(m_bp_T * mindelta) / chi_conductivity;
        } else if (fabs(eta_T) < SEPS && fabs(m_bp_T) > 0) { // m_bp_T is set
            eta_T = sq(m_bp_T * mindelta) / chi_conductivity;
        } else { // only eta is set
            m_bp_T = sqrt(eta_T * chi_conductivity) / mindelta;
        }
        fprintf(ferr, "Brinkman penalization params for the heat equation: eta_T=%g, m_bp_T=%g, minDelta=%g\n", eta_T, m_bp_T, mindelta);
    }else{
		eta_T = 1e+15;
		m_bp_T = 1e+8;
	}
}

// integration time step
double CFL_ARR = 0.3;//max changing of alpha in one timestep
event stability(i++){
  if (i % 100 == 0){
  	double dt_Arr_min = 1e+10;
  	foreach( reduction(min:dt_Arr_min)){
  	    double dt_arr = CFL_ARR / (KT(T[]) * dFR_dalpha(alpha_doc[]) + SEPS);
        if (dt_arr > dt_Arr_min) dt_Arr_min = dt_arr;
    }
  	dt = min(dt, dt_Arr_min);
  	fprintf(ferr, "dt_Arr=%g min(dt,dt_arr)=%g\n", dt_Arr_min, dt);
  }
}
