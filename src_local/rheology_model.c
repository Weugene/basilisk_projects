//The module is intended to solve the heat transfer equation, the polymerization effect and rheology changing
//it is coupled with navier-stokes/centered.h module
//
double Htr = 1;
//double Arrhenius_const = 4.53e+7;//1/s
//double Ea_by_R = 72900/8.314;// Kelvin
double Arrhenius_const = 10;//1/s
double Ea_by_R = 5;// Kelvin
double n_degree = 1.667;
double m_degree = 0.333;

//double Cp1 = 1100, Cp2 = 1006, Cp3 = 840;//J/(kg*K)
//double kappa1 = 1100, kappa2 = 0.02535, kappa3 = 1.11;//W/(m*K)
double Cp1 = 1, Cp2 = 1, Cp3 = 1;//J/(kg*K)
double kappa1 = 1, kappa2 = 0.01, kappa3 = 1;//W/(m*K)

void average(scalar A, const scalar f, const scalar fs, const double A1, const double A2, const double A3){
    foreach(){
        A[] = f[]*(A1 - A2) + A2 + fs[]*(A3 - A2);
    }
}
void advection_centered (scalar f, vector u, scalar df)
{
    foreach()
    df[] = ((f[] + f[-1,0])*u.x[] -
            (f[] + f[1,0])*u.x[1,0] +
            (f[] + f[0,-1])*u.y[] -
            (f[] + f[0,1])*u.y[0,1])/(2.*Delta);
    boundary ((scalar*) {df});
}

void advection_upwind (scalar f, vector u, scalar df)
{
    foreach()
    df[] = ((u.x[] < 0. ? f[] : f[-1,0])*u.x[] -
            (u.x[1,0] > 0. ? f[] : f[1,0])*u.x[1,0] +
            (u.y[] < 0. ? f[] : f[0,-1])*u.y[] -
            (u.y[0,1] > 0. ? f[] : f[0,1])*u.y[0,1])/Delta;
    boundary ((scalar*) {df});
}

event stability (i++) {
    dtmax = min(dtmax, 0.5*(h/L0)/(Arrhenius_const*exp(-Ea_by_R/TMAX)));
//    double cfl;
//    foreach_face(x, reduction (max:cfl)) {
//        cfl =
//    }
}

mgstats mgT;
scalar r[], thetav[];
scalar u_grad_scalar[], tmp[];
event end_timestep (i++,last){
//average(thetav, f, fs, rho1*Cp1, rho2*Cp2, rho3*Cp3);
foreach() thetav[] = f[]*(rho1*Cp1 - rho2*Cp2) + rho2*Cp2 + fs[]*(rho3*Cp3 - rho2*Cp2);
        foreach_face() kappa.x[] = f[]*(kappa1 - kappa2) + kappa2 + fs[]*(kappa3 - kappa2);
//    fprintf(stderr, "kappa, thetav done, grad will");

advection_centered(T, u, u_grad_scalar);
//    advection_upwind (T, u, u_grad_scalar);
//    fprintf(stderr, "grad adv done, tmp r will");
foreach() {
    tmp[] = Arrhenius_const*pow(1-alpha_doc[], n_degree)*exp(-Ea_by_R/T[]);
    r[] = Htr*rho1*f[]*tmp[] - thetav[]*u_grad_scalar[];
    //    printf("thetav=%g\n",thetav[]);
}
//    fprintf(stderr, "grad r tmp, diffusion will");
mgT = diffusion (T, dt, D = kappa, r = r, theta = thetav);

fprintf (stderr, "mg: i=%d t=%g p=%d u=%d T=%d\n", i, t, mgp.i, mgu.i, mgT.i); //number of iterations
advection_centered(alpha_doc, u, u_grad_scalar);
//    fprintf(stderr, "diffusion and advec done, alpha_doc evol will");
//    advection_upwind (alpha_doc, u, u_grad_scalar);
foreach() {
    //    alpha_doc[] += dt*(tmp[] - u_grad_scalar[]);
    alpha_doc[] = (alpha_doc[] + dt*(tmp[]*(1.0 + n_degree*alpha_doc[]/(1-alpha_doc[]+SEPS)) - u_grad_scalar[]))/(1 + dt*tmp[]*n_degree/(1-alpha_doc[]+SEPS));
    alpha_doc[] = clamp(alpha_doc[], 0.0, 1.0);
}
//    fprintf(stderr, "alpha_doc done");
boundary ((scalar*) {alpha_doc});
//    fprintf(stderr, "boundary...last");

}