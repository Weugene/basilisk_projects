/**
# Barotropic implicit time integration

This must be used in combination with the [multilayer solver](hydro.h).
It also includes Coriolis/linear friction terms i.e.
$$
\begin{aligned}
  \partial_t \left( h \mathbf{u} \right)_k + \mathbf{{\nabla}} \cdot \left(
  h \mathbf{u}  \mathbf{u} \right)_k & = - gh_k  \mathbf{{\nabla}} (\eta)
  {\color{blue}+ \mathbf{B}(h\mathbf{u})_k},
\end{aligned}
$$
with
$$
\mathbf{B} = \left(\begin{matrix}
-k & f \\
-f & -k
\end{matrix}\right)
$$
and $f$ the Coriolis parameter and $k$ the linear friction coefficient.

The discretisation scheme is inspired from [Popinet & Rickard,
2007](#popinet2007) (abbreviated as P & R below). 

We first remove the CFL constraint on barotropic waves and adjust the
defaults $\gamma_H$ and $\theta_H$ parameters. */

event defaults (i = 0)
{
  CFL_H = HUGE;
  gamma_H = theta_H = 0.5;
}

/**
## Coriolis and friction 

By default, the Coriolis and friction parameters ($f$ and $k$) are set
to zero. The user can change this by defining the macros `F0()` and
`K0()` before including the header file. */

#ifndef F0
# define F0() 0.
#endif
#ifndef K0
# define K0() 0.
#endif

/**
The Coriolis and friction terms are taken into account by linear
inversion of the $2 \times 2$ system
$$
(\mathbf{I} - \alpha_H\Delta t\mathbf{B})\hat{\mathbf{u}} = 
\left[\mathbf{I} + (1 - \alpha_H)\Delta t\mathbf{B}\right]\mathbf{u}^n -
2(1 - \gamma_H)\Delta t g\mathbf{\nabla}\eta^\star + \Delta t\mathbf{S}
$$
corresponding to equation (4) of P & R, where $\mathbf{S}$ includes
other source terns (e.g. advection) and $\eta^\star$ is $\eta^n$ in P
& R and $\eta^{n+1}$ here (as computed by explicit advection of $\eta$
[here](hydro.h#advection-and-diffusion)). */

double alpha_H = 0.5;

event coriolis (i++)
{
  foreach() {
    coord b0 = { - K0(), - K0() }, b1 = { F0(), -F0() };
    coord m0 = { 1. - alpha_H*dt*b0.x, 1. - alpha_H*dt*b0.y };
    coord m1 = { - alpha_H*dt*b1.x, - alpha_H*dt*b1.y };
    double det = m0.x*m0.y - m1.x*m1.y;
    foreach_layer()
      if (h[] > dry) {
        coord r, a;
	foreach_dimension() {
	  a.x = - (1. - gamma_H)*dt*(ha.x[] + ha.x[1])/(hf.x[] + hf.x[1] + dry);
	  r.x = u.x[] + (1. - alpha_H)*dt*(b0.x*u.x[] + b1.x*u.y[]) - 2.*a.x;
	}

	/**
	The intermediate velocity field is then obtained as
	$$
	\mathbf{u}^\star = \hat{\mathbf{u}} + 
	(1 - \gamma_H)\Delta tg\mathbf{\nabla}\eta^\star
	$$
	corresponding to equation (5) of P & R. */
	
	foreach_dimension()	  
	  u.x[] = (m0.y*r.x - m1.x*r.y)/det + a.x;
      }
  }
  boundary ((scalar *){u});  
}

/**
## Implicit free surface 

The free-surface elevation at the end of the timestep, $\eta^{n+1}$,
is computed by inversion of the Poisson--Helmholtz equation
$$
\eta^{n+1} - \theta_H\gamma_H\Delta t^2g\mathbf{\nabla}\cdot 
h^{n+1}\mathbf{\nabla}\eta^{n+1} = \eta^\star - \Delta t\mathbf{\nabla}\cdot
h^{n+1}\left[(1 - \theta_H)\mathbf{u}^n + \theta_H\mathbf{u}^\star\right]
$$
corresponding to equation (6) of P & R (see also eq. (58) of [Zijlema
& Stelling, 2005](#zijlema2005)).
*/

#include "poisson.h"

mgstats mgH;

event pressure (i++)
{

  /**
  We first compute
  $$
  \begin{aligned}
  (h\mathbf{u}^\star)_k & = h^{n+1}_k\left((1 - \theta_H)\mathbf{u}^n + 
                            \theta_H\mathbf{u}^\star\right)_k \\
  \alpha & = \sum h_k
  \end{aligned}
  $$
  */
  
  face vector alpha[];
  foreach_face() {
    alpha.x[] = 0.;
    foreach_layer() {
      double hl = h[-1] > dry ? h[-1] : 0.;
      double hr = h[] > dry ? h[] : 0.;
      double uf = hl > 0. || hr > 0. ? (hl*u.x[-1] + hr*u.x[])/(hl + hr) : 0.;
      hu.x[] += theta_H*hf.x[]*uf;
      alpha.x[] += hf.x[];
    }
    alpha.x[] *= 2.*fm.x[]/(cm[] + cm[-1]);
  }
  boundary ((scalar *){hu, alpha});
  
  scalar rhs[], etap[];
  for (int i = 0; i < nboundary; i++) {
    etap.boundary[i] = eta.boundary[i];
    etap.boundary_homogeneous[i] = eta.boundary_homogeneous[i];
  }

  /**
  The r.h.s. of (6) is computed as
  $$
  \frac{\eta^\star - \Delta t\mathbf{\nabla}\cdot\sum(h\mathbf{u}^\star)_k}
  {\lambda}
  $$
  with
  $$
  \lambda = -\theta_H\gamma_H \Delta t^2 g
  $$
  */
  
  double H1 = 0., v1 = 0.;
  foreach(reduction(+:H1) reduction(+:v1)) {
    rhs[] = etap[] = eta[];
    foreach_layer() {
      H1 += dv()*h[];
      foreach_dimension()
	rhs[] -= dt*(hu.x[1] - hu.x[])/Delta;
    }
    rhs[] /= - theta_H*gamma_H*G*sq(dt);
    v1 += dv();
  }
  boundary ({etap});

  /**
  The Poisson--Helmholtz equation is solved as
  $$
  \lambda\eta^{n+1} + \mathbf{\nabla}\cdot\alpha\mathbf{\nabla}\eta^{n+1} 
  = \text{rhs}
  $$
  */
  
  const scalar lambda[] = - 1./(theta_H*gamma_H*G*sq(dt));
  mgH = poisson (etap, rhs, alpha, lambda = lambda,
		 tolerance = TOLERANCE*H1/(v1*G*sq(dt)));

  /**
  And finally the barotropic term is added to the face-weighted
  accelerations and fluxes, $ha_{i+1/2}^{n+1}$ and
  $hu_{i+1/2}^{n+1}$ (step 4. of P & R, page 7).*/

  foreach_face() {
    double ax = - fm.x[]*G*(etap[] - etap[-1])/((cm[] + cm[-1])*Delta/2.);
    foreach_layer() {
      ha.x[] = hf.x[]*ax;
      hu.x[] += theta_H*gamma_H*dt*ha.x[];
    }
  }
  boundary ((scalar *){ha, hu});
}

/**
## References

~~~bib
@Article{popinet2007,
  author = 	 {S. Popinet and G. Rickard},
  title = 	 {A tree-based solver for adaptive ocean modelling},
  journal = 	 {Ocean Modelling},
  year = 	 {2007},
  number =       {16},
  pages =        {224-249},
  url =          {http://gfs.sf.net/ocean.pdf}
}

@article{zijlema2005,
  title={Further experiences with computing non-hydrostatic free-surface flows involving water waves},
  author={Zijlema, Marcel and Stelling, Guus S},
  journal={International Journal for Numerical Methods in Fluids},
  volume={48},
  number={2},
  pages={169--197},
  year={2005},
  publisher={Wiley Online Library}
}
~~~
*/
