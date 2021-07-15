/**
# Three-phase interfacial flows

This file helps setup simulations for flows of two fluids separated by
an interface (i.e. immiscible fluids) with solid obstacles. It is typically used in
combination with a [Navier--Stokes solver](navier-stokes/centered.h).

The interface between the fluids is tracked with a Volume-Of-Fluid
method. The volume fraction $f$ and $fs$ leads to next interpretation of averaged value of A =(\rho, \mu, \kappa)
\begin{table}[]
\begin{tabular}{lll}
f & fs & A    \\
1 & 1  & A\_3 \\
0 & 1  & A\_3 \\
1 & 0  & A\_1 \\
0 & 0  & A\_2
\end{tabular}
\end{table}
 The above definition of variables leads to a specific definition of the averaged value A.
The densities and dynamic viscosities for fluid 1 and 2 are *rho1*,
*mu1*, *rho2*, *mu2*, respectively, in a solid rho3 and mu3.*/

#include "vof.h"
double VOF_cutoff = 0.01;
scalar f[], fs[], * interfaces = {f}, * interfaces_all = {f,fs};
double mu0 = 0, rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0., rho3 = 1., mu3 = 0.;
double kappa1 = 0, kappa2 = 0, kappa3 = 0;//W/(m*K)
double Cp1 = 0, Cp2 = 0, Cp3 = 0;
int N_smooth = 3;
/**
Auxilliary fields are necessary to define the (variable) specific
volume $\alpha=1/\rho$ as well as the cell-centered density. */
face vector alphav[];
scalar rhov[];
#ifdef HEAT_TRANSFER
	face vector kappav[];
#endif

event defaults (i = 0) {
    alpha = alphav;
    rho = rhov;

    /**
    If the viscosity and conductivity are non-zero, we need to allocate the face-centered
    viscosity and conductivity fields. */

    if (mu1 || mu2) //?
        mu = new face vector;
#ifdef HEAT_TRANSFER
    if (kappa1 || kappa2) //?
        kappa = new face vector;
#endif
    /**
    We add the interface to the default display. */

    display ("draw_vof (c = 'f');");
}

/**
The density, viscosity and conductivity are defined using arithmetic averages by
default. The user can overload these definitions to use other types of
averages (i.e. harmonic).
Usually, it is assumed that mu1 is variable, mu2 and mu3 are not. For simplisity mu3=mu2
 */
#ifndef var_hom
#define var_hom(f, fs, A1, A2, A3) ((1.0 - clamp(fs,0.,1.))*(A2 + (A1 - A2)*clamp(f,0.,1.)) + clamp(fs,0.,1.)*A3)
#endif

#ifndef var_harm
#define var_harm(f, fs, A1, A2, A3) (1.0/(  (1.0 - clamp(fs,0.,1.))*( clamp(f,0.,1.)*(1.0/(A1) - 1.0/(A2)) + 1.0/(A2)) + clamp(fs,0.,1.)/(A3)))

#endif

#ifndef rho
#define rho(f, fs) var_hom(f, fs, rho1, rho2, rho3)
//#define rho(f, fs) ((1.0 - clamp(fs,0.,1.))*(rho2 + (rho1 - rho2)*clamp(f,0.,1.)) + clamp(fs,0.,1.)*rho3)
#endif

#ifndef kappav
//    #define kappav(f, fs) var_harm(f, fs, kappa1, kappa2, kappa3)
    #define kappav(f, fs) ((1.0 - clamp(fs,0.,1.))*(kappa2 + (kappa1 - kappa2)*clamp(f,0.,1.)) + clamp(fs,0.,1.)*kappa3)
#endif
#if REACTION_MODEL != NO_REACTION_MODEL
    scalar alpha_doc[];
#endif
scalar T[];
/**
# Variable rheology models
 $$\mu = \mu_1 \exp(\frac{E_\eta}{RT})(\frac{\alpha_{gel}}{\alpha_{gel}-\alpha})^f(\alpha, T)$$
 $$f(\alpha, T) = A + B \alpha$$

 $$\mu = \mu_1 \exp(\frac{E_\eta}{RT}+\chi \alpha^2)$$
 **/
double Eeta_by_Rg = 0.1; //Kelvin
double chi = 1;
double alpha_gel = 0.8;
double mu_eff = 0;
#ifndef fpol
    #define fpol(alpha_doc, T) A*alpha_doc + B
#endif


#ifndef muf1
    #if REACTION_MODEL != NO_REACTION_MODEL
        //#define mupol(alpha_doc, T) (mu0*exp(Eeta_by_Rg/(T))*pow(alpha_gel/(alpha_gel-alpha_doc), fpol(alpha_doc, T)))
        #define muf1(alpha_doc, T) (mu0*exp(Eeta_by_Rg/T + chi*alpha_doc))
    #else
        #define muf1(alpha_doc, T) (mu0*exp(Eeta_by_Rg/T))
    #endif
#endif

#ifndef mu
//#define mu(f, fs, alpha_doc, T) (mu2 + (mu1 - mu2)*clamp(f,0.,1.))
//#define mu(f, fs, alpha_doc, T) ((1.0 - clamp(fs,0.,1.))*(mu2 + (muf1(alpha_doc, T) - mu2)*clamp(f,0.,1.)) + clamp(fs,0.,1.)*mu3)
#define mu(f, fs, alpha_doc, T) (1.0/(  (1.0 - clamp(fs,0.,1.))*( clamp(f,0.,1.)*(1.0/muf1(alpha_doc, T) - 1.0/mu2) + 1.0/mu2) + clamp(fs,0.,1.)/mu3))
#endif

/**
We have the option of using some "smearing" of the density/viscosity
jump. */

#ifdef FILTERED
    scalar sf1[], sf2[];
#else
    #define sf1 f
    #define sf2 fs
#endif
void filter_scalar(scalar f, scalar sf){
#if dimension <= 2
    foreach()
        sf[] = (4.*f[] +
             2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
             f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
#else // dimension == 3
    foreach()
        sf[] = (8.*f[] +
            4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
            2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] +
                f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
                f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
            f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
            f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
#endif
#if TREE
    sf.prolongation = refine_bilinear;
    boundary ({sf});
#endif
}

event tracer_advection (i++) {
	scalar sf_s[];
    /**
    When using smearing of the density jump, we initialise *sf* with the
    vertex-average of *f*. */

#ifndef sf1
		filter_scalar(f, sf1);
		for (int i_smooth=2; i_smooth<=N_smooth; i_smooth++){
			filter_scalar(sf1, sf_s);
			foreach() sf1[] = sf_s[];
			boundary ({sf1});
		}
#endif

#ifndef sf2
    filter_scalar(fs, sf2);
		for (int i_smooth=2; i_smooth<=N_smooth; i_smooth++){
			filter_scalar(sf2, sf_s);
			foreach() sf2[] = sf_s[];
			boundary ({sf2});
		}
#endif

}

#define n_mu_eff 4
#define norm_mu_eff (sq(n_mu_eff) + 3*n_mu_eff + 2)
event properties (i++) {
#ifdef DAMP_CAPILLARY_WAVE
    scalar f_very_smooth[], sf_s[];
    filter_scalar(f, f_very_smooth);
		for (int i_smooth=2; i_smooth<=5; i_smooth++){
			filter_scalar(f_very_smooth, sf_s);
			foreach() f_very_smooth[] = sf_s[];
			boundary ({f_very_smooth});
		}
#endif
    foreach_face() {
        double ff1 = (sf1[] + sf1[-1])/2.; //liquid
        double ff2 = (sf2[] + sf2[-1])/2.; //solid
        alphav.x[] = fm.x[]/rho(ff1, ff2);
        if (mu1 || mu2) {
            face vector muv = mu;
            double Tf = (T[] + T[-1])/2.;
#if REACTION_MODEL != NO_REACTION_MODEL
            double alpha_doc_f = (alpha_doc[] + alpha_doc[-1])/2.;
            muv.x[] = fm.x[]*mu(ff1, ff2, alpha_doc_f, Tf);
#else
            muv.x[] = fm.x[]*mu(ff1, ff2, 0, Tf);
#endif
#ifdef DAMP_CAPILLARY_WAVE
					 muv.x[] += fm.x[]*mu_eff*(ff1) * pow(1 - ff1, n_mu_eff)*norm_mu_eff;
#endif
        }
#ifdef HEAT_TRANSFER
        if (kappa1 || kappa2) {
            face vector kappav = kappa;
            kappav.x[] = fm.x[]*kappav(ff1, ff2);
        }
#endif
    }
    foreach()
        rhov[] = cm[]*rho(sf1[], sf2[]); //? alphav.x and rhov are not consistent - All do so
#if TREE
    sf1.prolongation = fraction_refine; //after changing we restore prolongation operator
    sf2.prolongation = fraction_refine;
    boundary ({sf1,sf2});
#endif
}

double give_etas(double m_bp, double mindelta, double nu_bp){
	return sq(m_bp * mindelta) / nu_bp;
}
double give_mbp(double eta_s, double mindelta, double nu_bp){
	return sqrt(eta_s * nu_bp) / mindelta;
}
event properties (i < 10){
    double nu_bp = mu1/rho1, mindelta=1e+10;
    foreach( reduction(min:mindelta) ){
        if (Delta < mindelta) mindelta = Delta;
    }

    if (nu_bp) {
        if (fabs(eta_s) > SEPS && fabs(m_bp) > 0) { // m_bp has higher priority
            eta_s = sq(m_bp * mindelta) / nu_bp;
        } else if (fabs(eta_s) < SEPS && fabs(m_bp) == 0) { // nothing is set
            m_bp = 1;
            eta_s = sq(m_bp * mindelta) / nu_bp;
        } else if (fabs(eta_s) < SEPS && fabs(m_bp) > 0) { // m_bp is set
            eta_s = sq(m_bp * mindelta) / nu_bp;
        } else { // only eta is set
            m_bp = sqrt(eta_s * nu_bp) / mindelta;
        }
        fprintf(ferr, "Brinkman penalization params: eta_s=%g, m_bp=%g, minDelta=%g\n", eta_s, m_bp, mindelta);
    }
}
