/**
# Surface tension

Surface tension can be expressed as the interfacial force density
$$
\phi\nabla f
$$
with $f$ the volume fraction describing the interface and the potential
$$
\phi = \sigma\kappa
$$
with $\sigma$ the (constant) surface tension coefficient and $\kappa$
the interface mean curvature. */
#include "curvature.h"

attribute {
    scalar phi;
}

event defaults (i = 0) {
    if (is_constant(a.x)) {
        a = new face vector;
        foreach_face()
        a.x[] = 0.;
        boundary ((scalar *){a});
    }
}

////////////////////////////////////////////////////////////////////



/**
The surface tension coefficient is associated to each VOF tracer. */

double sigma12 = 0., sigma23 = 0., sigma13 = 0.;
double sigma1 = 0., sigma2 = 0., sigma3 = 0.;

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

/**
## Definition of the potential

We overload the acceleration event to define the potential
$\phi=\sigma\kappa$.<br/>

### This is where we make the change
We first calculate the curvature $\kappa$, instead of $sigma\kappa$. Then we change this field
based on the interface location.
*/

event acceleration (i++)
{
    /**
    We check for all existing interfaces_all (even which does not move) for which $\sigma$ is non-zero. */
    for (scalar f in interfaces_all){
        /**
        If $\phi$ is already allocated, an error is thrown. (Fix me!!)
        */
        scalar phi = f.phi;
        if (phi.i){
            fprintf(ferr, "Phi already exists! The current ThreePhase model cannot handle that.\n");
        }
        phi = new scalar;
        curvature (f, phi, 1.0, add = false);
        f.phi = phi;
    }
    /**
    Right now, phi only means $\kappa$. We need to change it to \sigma\kappa based on the interface locations.<br/>
    We refer back to the Schematic:
    <figure>
    <p align="center">
      <img src="https://www.dropbox.com/s/3l0rvxwh6cgrt7g/SchematicBBasilisk.png?dl=1" width="40%">
      <figcaption><p align="center">Schematic: Details of the VOF fields and interfaces.</figcaption>
    </figure>
    * If both phi1[] and phi2[] are defined, then it is the interface between phase 2 and 3.
    * If only phi1[] is defined (but phi2[]) is not, then it is the interface between 1 and 3.
    * If only phi2[] is defined (but phi1[]) is not, then it is the interface between 1 and 2.
    We modify $phi$ accordingly.
    */
    scalar phi1 = f.phi, phi2 = fs.phi;
    foreach(){ // To do: Use a better if-else look. Also check if for (scalar f in interfaces) can be used somehow.
        if(phi1[] < nodata && phi2[] < nodata){
            phi1[] *= sigma23/2.;
            phi2[] *= sigma23/2.;
        } else {
            if(phi1[] < nodata){
                phi1[] *= sigma13;
            } else {
                if (phi2[] < nodata){
                    phi2[] *= sigma12;
                }
            }
        }
    }
///////////////////////////////////////////////////////
    /**
    We check for all existing interfaces for which $\phi$ is allocated. The
    corresponding volume fraction fields will be stored in *list*. */

    scalar * list = NULL;
    for (scalar f in interfaces_all)
        if (f.phi.i) {
            list = list_add (list, f);

            /**
            To avoid undeterminations due to round-off errors, we remove
            values of the volume fraction larger than one or smaller than
            zero. */

            foreach()	f[] = clamp (f[], 0., 1.);
            boundary ({f});
        }

    /**
    On trees we need to make sure that the volume fraction gradient
    is computed exactly like the pressure gradient. This is necessary to
    ensure well-balancing of the pressure gradient and interfacial force
    term. To do so, we apply the same prolongation to the volume
    fraction field as applied to the pressure field. */

#if TREE
    for (scalar f in list)
    f.prolongation = p.prolongation;
  boundary (list);
#endif

    /**
    Finally, for each interface for which $\phi$ is allocated, we
    compute the interfacial force acceleration
    $$
    \phi\mathbf{n}\delta_s/\rho \approx \alpha\phi\nabla f
    $$
    */

    face vector ia = a;
    foreach_face()
    for (scalar f in list)
        if (f[] != f[-1]) {

            /**
            We need to compute the potential *phif* on the face, using its
            values at the center of the cell. If both potentials are
            defined, we take the average, otherwise we take a single
            value. If all fails we set the potential to zero: this should
            happen only because of very pathological cases e.g. weird
            boundary conditions for the volume fraction. */

            scalar phi = f.phi;
            double phif =
                    (phi[] < nodata && phi[-1] < nodata) ?
                    (phi[] + phi[-1])/2. :
                    phi[] < nodata ? phi[] :
                    phi[-1] < nodata ? phi[-1] :
                    0.;

            ia.x[] += alpha.x[]/fm.x[]*phif*(f[] - f[-1])/Delta;
        }

    /**
    On trees, we need to restore the prolongation values for the
    volume fraction field. */

#if TREE
    for (scalar f in list)
    f.prolongation = fraction_refine;
  boundary (list);
#endif

    /**
    Finally we free the potential fields and the list of volume
    fractions. */

    for (scalar f in list) {
        scalar phi = f.phi;
        delete ({phi});
        f.phi.i = 0;
    }
    free (list);
}
