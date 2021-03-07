//The module is intended to solve the heat transfer equation, the polymerization effect and rheology changing
//it is coupled with navier-stokes/centered.h module
//
#define HEAT_TRANSFER
#include "three-phase-rheology.h"
#include "diffusion-weugene.h"
//const scalar const_temp_solid[] = 1.1;
(const) scalar T_target = unity;
double eta_T = 0;
int m_bp_T = 0;
bool stokes_heat = false;
double chi_conductivity = 0;
double Htr = 1;
double Arrhenius_const = 10;//1/s
double Ea_by_R = 5;// Kelvin
double n_degree = 1.667;
double m_degree = 0.333;

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
#define REACTION_MODEL_NON_AUTOCATALYTIC 1
#define REACTION_MODEL_N_ORDER_AUTOCATALYTIC 2
#define REACTION_MODEL_PROUT_TOMPKINS_AUTOCATALYTIC 3

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
#else
#define FR(alpha_doc) ( pow(1 - alpha_doc, n_degree) )
#define dFR_dalpha(alpha_doc) ( -n_degree * pow(1 - alpha_doc, n_degree - 1) )
#define GENERAL_METHOD 0
#endif

mgstats mgT;
event end_timestep (i++){
    scalar r[], thetav[], tmp[], beta[];
    foreach() thetav[] =  var_hom(f[], fs[], rho1*Cp1, rho2*Cp2, rho3*Cp3);
    foreach_face() kappav.x[] = kappav(f[], fs[]);
    boundary ({kappav, thetav});
    //The advection step
    if (!stokes_heat) {
        advection((scalar *) {T}, uf, dt);
    }
    //source and conductivity terms. f[] * (1 - fs[]) multiplications means gas and solids can't produce heat
    //solids play role in the conduction process

    if (fabs(Htr) > SEPS){
//        fprintf(ferr, "+++ %d %g %g\n", GENERAL_METHOD, Htr, SEPS);

        foreach() {
            #if GENERAL_METHOD == 0
                tmp[] = f[] * (1 - fs[]) * Arrhenius_const * exp(-Ea_by_R / T[]) * pow(1 - alpha_doc[], n_degree);
                r[] = Htr * rho1 * tmp[] * (1.0 - Eeta_by_Rg/T[]);
                beta[] = Htr * rho1 * tmp[] * Eeta_by_Rg/(T[]*T[]);
            #else  //General case
                tmp[] = Htr * rho1 * f[] * (1 - fs[]) * (1 - fs[]) * FR(alpha_doc[]);
                r[] = tmp[] * ( KT(T[]) - dKT_dT(T[]) * T[]);
                beta[] = tmp[] * dKT_dT(T[]);
            #endif
                //Penalization terms are:
                r[] += fs[] * thetav[] * T_target[] / eta_T;
                beta[] += -fs[] * thetav[] / eta_T;
        }
    }else {
        foreach() {
            r[] = fs[] * thetav[] * T_target[] / eta_T;
            beta[] = -fs[] * thetav[] / eta_T;
            tmp[] = f[] * (1 - fs[]) * KT(T[]) * FR(alpha_doc[]);
        }
    }
    boundary((scalar *){r, beta, tmp});
//        fprintf(stderr, "diffusion");
    if (constant(kappa.x) != 0.) {
        mgT = diffusion(T, dt, D = kappa, r = r, beta = beta, theta = thetav);
        fprintf (stderr, "mgT: i=%d t=%g dt=%g num of iterations T=%d\n", i, t, dt, mgT.i); //number of iterations
    }else{
        foreach(){
            T[] = (thetav[] * T[] + dt * r[])/(thetav[] - dt * beta[]);
        }
    }
    boundary((scalar *){T});
    // advection term of kinetic equation is solved implicitly
    // due to a linearized source term
    advection ((scalar *){alpha_doc}, uf, dt);
    foreach() {
        #if GENERAL_METHOD == 777
			alpha_doc[] = 1.0 - pow(pow(fabs(1 - alpha_doc[]), 1 - n_degree) + (n_degree - 1.0) * dt * f[] * (1 - fs[]) * KT(T[]), 1.0 - n_degree);//direct integration from t to t + dt at fixed T
            //alpha_doc[] = (alpha_doc[] + dt * (tmp[] * (1.0 + n_degree*alpha_doc[] / (1 - alpha_doc[] + SEPS))))/(1 + dt * tmp[] * n_degree / (1 - alpha_doc[] + SEPS));//numerical solution
        #else
            alpha_doc[] = (alpha_doc[] + dt * f[] * KT(T[]) * (FR(alpha_doc[]) - dFR_dalpha(alpha_doc[]) * alpha_doc[])) /
                          (1 - dt * f[] * KT(T[]) * dFR_dalpha(alpha_doc[]));
        #endif
        alpha_doc[] = clamp(alpha_doc[], 0.0, 1.0);
    }
    boundary ((scalar*) {alpha_doc});
}

// In first 10 steps, eta_T and m_bp_T will be adjusted
event properties (i += 1; i < 10){
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
        fprintf(ferr, "Brinkman penalization params for the heat equation: eta_T=%g, m_bp_T=%d, minDelta=%g\n", eta_T, m_bp_T, mindelta);
    }else{
		eta_T = 1e+15;
		m_bp_T = (int) 1e+8;
	}
}

// integration time step
double CFL_ARR = 0.01;//max changing of alpha in one timestep
event stability(i++){
	double maxT = -1e+10, mindelta = 0;
	foreach( reduction(min:mindelta) reduction(max:maxT)){
		if (T[] > maxT) maxT = T[];
    }
	double dt_Arr = CFL_ARR /(Arrhenius_const * exp(-Ea_by_R / maxT));
	dt = min(dt, dt_Arr);
	fprintf(ferr, "dt_Arr=%g dt_cur=%g", dt_Arr, dt);
}
