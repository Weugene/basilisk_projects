#include "curvature.h"

attribute {
    scalar phi;
    double sigma;
}

scalar phi1[];
scalar phi2[];
double sigma_sum = 0;
double sigma12 = 0., sigma23 = 0., sigma13 = 0.;
double sigma3 = 0.;

event defaults (i = 0) {
    if (is_constant(a.x)) {
        a = new face vector;
        foreach_face() a.x[] = 0.;
        boundary ((scalar *){a});
    }
    sigma_sum = sigma12 + sigma13 + sigma23;
    f.phi = phi1;
    fs.phi = phi2;
//    fprintf(stderr, "sigma12=%g sigma13=%g sigma23=%g", sigma12, sigma13, sigma23);
    if (sigma_sum) {
        f.sigma = 0.5*(sigma13 + sigma12 - sigma23);
        fs.sigma = 0.5*(- sigma13 + sigma12 + sigma23);
        sigma3 = 0.5*(sigma13 - sigma12 + sigma23);
        double costheta = (sigma23 - sigma12)/sigma13;
        fprintf(stderr, "cos theta=%g theta=%g\n", costheta, 180*acos(costheta)/3.141592);
        fprintf(stderr, "sigma1=%g sigma2=%g sigma3=%g\n", f.sigma, fs.sigma, sigma3);
    }
}

/**
## Stability condition

The surface tension scheme is time-explicit so the maximum timestep is
the oscillation period of the smallest capillary wave.
$$
T = \sqrt{\frac{\rho_{m}\Delta_{min}^3}{\pi\sigma}}
$$
with $\rho_m=(\rho_1+\rho_2)/2.$ and $\rho_1$, $\rho_2$ the densities
on either side of the interface. */

event stability (i++) {
    if (!sigma_sum) return 0; //???
    /**
    We first compute the minimum and maximum values of $\alpha/f_m =
    1/\rho$, as well as $\Delta_{min}$. */

    double amin = HUGE, amax = -HUGE, dmin = HUGE;
    foreach_face (reduction(min:amin) reduction(max:amax) reduction(min:dmin)) {
        if (alpha.x[]/fm.x[] > amax) amax = alpha.x[]/fm.x[];
        if (alpha.x[]/fm.x[] < amin) amin = alpha.x[]/fm.x[];
        if (Delta < dmin) dmin = Delta;
    }
    double rhom = (1./amin + 1./amax)/2.;
    /**
    The maximum timestep is set using the sum of surface tension
    coefficients. */

    double sigma = sigma13+sigma23+sigma12;

    if (sigma) {
        double dt = sqrt (rhom*cube(dmin)/(pi*sigma));
        if (dt < dtmax)
            dtmax = dt;
    }
}

event acceleration (i++)
{
    if (!sigma_sum) return 0; //???
    /**
    We check for all existing interfaces_all (even which does not move) for which $\sigma$ is non-zero. */

//    for (scalar f in interfaces_all){
//        scalar phi = f.phi;
//        if (phi1.i && phi2.i){
//            fprintf(ferr, "Phi already exists! The current ThreePhase model cannot handle that.\n");
//        }
//        phi = new scalar;
//        curvature (f, phi, 1.0, add = false);
//        f.phi = phi;
////        foreach(){
////            fprintf(stderr, "kappa[%s] = %g", f.name, phi[]);
////        }
//    }




    curvature (f, phi1, 1.0, add = false);
    curvature (fs, phi2, 1.0, add = false);

//        foreach(){
//            fprintf(stderr, "kappa[%s] = %g", f.name, phi[]);
//        }




    /**
    On trees we need to make sure that the volume fraction gradient
    is computed exactly like the pressure gradient. This is necessary to
    ensure well-balancing of the pressure gradient and interfacial force
    term. To do so, we apply the same prolongation to the volume
    fraction field as applied to the pressure field. */

#if TREE
    for (scalar f in interfaces_all) //????
    f.prolongation = p.prolongation;
  boundary (interfaces_all);
#endif

    /**
    Finally, for each interface for which $\phi$ is allocated, we
    compute the interfacial force acceleration
    $$
    \phi\mathbf{n}\delta_s/\rho \approx \alpha\phi\nabla f
    $$
    */

//    face vector ia = a;
//    phi1 = f.phi;
//    phi2 = fs.phi;
//    foreach_face() {
//        for (scalar f in interfaces_all) {
//            scalar phi = f.phi;
//            ia.x[] += (phi[] < nodata) ? alpha.x[] / fm.x[] * f.sigma * phi[] * (f[1] - f[-1])/(2*Delta) : 0;
//        }
//        ia.x[] += (phi1[] < nodata && phi2[] < nodata) ? alpha.x[] / fm.x[] * sigma3 * (phi1[] + phi2[]) * (f[1] - f[-1] + fs[1] - fs[-1]) /(2*Delta) : 0; //gas
//    }
    face vector ia = a;
    foreach_face() {
        for (scalar f in interfaces_all) {
            scalar phi = f.phi;
            ia.x[] += (phi[] < nodata) ? alpha.x[] / fm.x[] * f.sigma * phi[] * (f[1] - f[-1])/(2*Delta) : 0;
//            if (phi[] && phi[] <nodata) fprintf(stderr, "ia=%g phi=%g f.sigma=%g\n", ia.x[], phi[], f.sigma);
        }
        ia.x[] += (phi1[] < nodata && phi2[] < nodata) ? alpha.x[] / fm.x[] * sigma3 * (phi1[] + phi2[]) * (f[1] - f[-1] + fs[1] - fs[-1]) /(2*Delta) : 0; //gas
//
    }
    /**
    On trees, we need to restore the prolongation values for the
    volume fraction field. */

#if TREE
    for (scalar f in interfaces_all) //?????
        f.prolongation = fraction_refine;
  boundary (interfaces_all);
#endif

}
