#ifndef BASILISK_HEADER_9
#define BASILISK_HEADER_9
#line 1 "./../src_local/rheology_model.h"
//The module is intended to solve the heat transfer equation, the polymerization effect and rheology changing
//it is coupled with navier-stokes/centered.h module
//
#include "../src_local/three-phase-rheology.h"
#include "diffusion.h"
scalar u_grad_scalar[];
double Htr = 1;
double Arrhenius_const = 10;//1/s
double Ea_by_R = 5;// Kelvin
double n_degree = 1.667;
double m_degree = 0.333;
double Cp1 = 1, Cp2 = 1, Cp3 = 1;//J/(kg*K)

#undef SEPS
#define SEPS 1e-10

//#ifndef AverA
//#define AverA(f, fs, A1, A2, A3)  (clamp((f-fs),0.,1.)*(A1 - A2) + A2 + clamp(fs,0.,1.)*(A3 - A2))
////#define AverA(f, fs, A1, A2, A3)  ( 1./(  clamp((f-fs),0.,1.)*( 1./A1 - 1./A2 ) + 1./A2 + clamp(fs,0,1)*(1./A3 - 1./A2)  ))
//#endif
//void average(scalar A, const scalar f, const scalar fs, const double A1, const double A2, const double A3){
//    foreach(){
//        A[] = f[]*(A1 - A2) + A2 + fs[]*(A3 - A2);
//    }
//}
void advection_centered (scalar f, vector u, scalar df)
{
    foreach()
    df[] = ((f[] + f[-1,0])*u.x[] -
            (f[] + f[1,0])*u.x[1,0] +
            (f[] + f[0,-1])*u.y[] -
            (f[] + f[0,1])*u.y[0,1])/(2.*Delta);
    boundary ((scalar*) {df});
}

//void advection_upwind (scalar f, vector u, scalar df)
//{
//    foreach()
//    df[] = ((u.x[] < 0. ? f[] : f[-1,0])*u.x[] -
//            (u.x[1,0] > 0. ? f[] : f[1,0])*u.x[1,0] +
//            (u.y[] < 0. ? f[] : f[0,-1])*u.y[] -
//            (u.y[0,1] > 0. ? f[] : f[0,1])*u.y[0,1])/Delta;
//    boundary ((scalar*) {df});
//}

event stability (i++) {
//    dtmax = min(dtmax, 0.5*(h/L0)/(Arrhenius_const*exp(-Ea_by_R/TMAX)));
//    double cfl;
//    foreach_face(x, reduction (max:cfl)) {
//        cfl =
//    }
}

mgstats mgT;
event end_timestep (i++){


    scalar r[], thetav[], tmp[], beta[];
    foreach() thetav[] = f[]*(rho1*Cp1 - rho2*Cp2) + rho2*Cp2 + fs[]*(rho3*Cp3 - rho2*Cp2);
//        foreach_face() kappav.x[] = kappav(f[], fs[]);
//        boundary ({kappav, thetav});
    advection_centered(T, u, u_grad_scalar);
    //    advection_upwind (T, u, u_grad_scalar);
//        fprintf(stderr, " tmp r ");
    foreach()
    {
        tmp[] = Arrhenius_const * pow(1 - alpha_doc[], n_degree) * exp(-Ea_by_R / T[]);
        r[] = Htr * rho1 * f[] * tmp[]*(1.0 - Eeta_by_Rg/T[]) - thetav[] * u_grad_scalar[];
        beta[] = Htr * rho1 * f[] * tmp[]*Eeta_by_Rg/(T[]*T[]);
    }
//        fprintf(stderr, "diffusion");
    if (constant(kappa.x) != 0.) {
        mgT = diffusion(T, dt, D = kappa, r = r, beta = beta, theta = thetav);
        fprintf (stderr, "mg: i=%d t=%g dt=%g p=%d u=%d T=%d\n", i, t, dt, mgp.i, mgu.i, mgT.i); //number of iterations
    }else{
        foreach(){
            T[] += dt*r[];
        }
    }

    advection_centered(alpha_doc, u, u_grad_scalar);
//        advection_upwind (alpha_doc, u, u_grad_scalar);
    foreach() {
        //    alpha_doc[] += dt*(tmp[] - u_grad_scalar[]);
        alpha_doc[] = (alpha_doc[] + dt*(tmp[]*(1.0 + n_degree*alpha_doc[]/(1-alpha_doc[]+SEPS)) - u_grad_scalar[]))/(1 + dt*tmp[]*n_degree/(1-alpha_doc[]+SEPS));
        alpha_doc[] = clamp(alpha_doc[], 0.0, 1.0);
    }

    boundary ((scalar*) {alpha_doc});
}
#endif
