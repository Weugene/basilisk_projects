#ifndef BASILISK_HEADER_6
#define BASILISK_HEADER_6
#line 1 "./../src_local/centered-weugene.h"
/**
# Incompressible Navier--Stokes solver (centered formulation)

We wish to approximate numerically the incompressible,
variable-density Navier--Stokes equations
$$
\partial_t\mathbf{u}+\nabla\cdot(\mathbf{u}\otimes\mathbf{u}) = 
\frac{1}{\rho}\left[-\nabla p + \nabla\cdot(2\mu\mathbf{D})\right] + 
\mathbf{a}
$$
$$
\nabla\cdot\mathbf{u} = 0
$$
with the deformation tensor 
$\mathbf{D}=[\nabla\mathbf{u} + (\nabla\mathbf{u})^T]/2$.

The scheme implemented here is close to that used in Gerris ([Popinet,
2003](/src/references.bib#popinet2003), [Popinet,
2009](/src/references.bib#popinet2009), [Lagrée et al,
2011](/src/references.bib#lagree2011)).

We will use the generic time loop, a CFL-limited timestep, the
Bell-Collela-Glaz advection scheme and the implicit viscosity
solver.
There are 2 ways of solid treatment: (i) embedded boundaries method (EBM), (ii) Brinkman Penalization method (BPM).
If embedded boundaries are used, a different scheme is used for viscosity.
The advantages of EBM are high accuracy and MPI parallelization, however, it works only for \textbf{one phase flow}.
If Brinkman Penalization method for solid treatment is used,
then a penalization term is added implicitly with the viscous term.
In this case $\mathbf{a} = -\frac{\chi}{\eta_s}\left(\mathbf{u} - \mathbf{U_t} \right)$
BPM can work for multiphase flow with MPI parallelization with sufficient accuracy.
Detailed information about comparison you can find here.
 */

#include "run.h"
#include "timestep.h"
#include "bcg.h"
#if EMBED
#include "viscosity-embed.h"
#else
#ifndef BRINKMAN_PENALIZATION
#include "viscosity.h"
#else
#include "viscosity-weugene.h"
#include "utils-weugene.h"
#endif
#endif

/**
The primary variables are the centered pressure field $p$ and the
centered velocity field $\mathbf{u}$. The centered vector field
$\mathbf{g}$ will contain pressure gradients and acceleration terms.

We will also need an auxilliary face velocity field $\mathbf{u}_f$ and
the associated centered pressure field $p_f$. */

scalar p[];
vector u[], g[];
scalar pf[];
face vector uf[];

/**
In the case of variable density, the user will need to define both the
face and centered specific volume fields ($\alpha$ and $\alpha_c$
respectively) i.e. $1/\rho$. If not specified by the user, these
fields are set to one i.e. the density is unity.

Viscosity is set by defining the face dynamic viscosity $\mu$; default
is zero.

The face field $\mathbf{a}$ defines the acceleration term; default is
zero.

The statistics for the (multigrid) solution of the pressure Poisson
problems and implicit viscosity are stored in *mgp*, *mgpf*, *mgu*
respectively. 

If *stokes* is set to *true*, the velocity advection term
$\nabla\cdot(\mathbf{u}\otimes\mathbf{u})$ is omitted. This is a
reference to [Stokes flows](http://en.wikipedia.org/wiki/Stokes_flow)
for which inertia is negligible compared to viscosity. */

(const) face vector mu = zerof, a = zerof, aold = zerof, alpha = unityf, kappa = zerof;
(const) scalar rho = unity;
mgstats mgp, mgpf, mgu;
bool stokes = false;

/**
## Boundary conditions

For the default symmetric boundary conditions, we need to ensure that
the normal component of the velocity is zero after projection. This
means that, at the boundary, the acceleration $\mathbf{a}$ must be
balanced by the pressure gradient. Taking care of boundary orientation
and staggering of $\mathbf{a}$, this can be written */

#if EMBED
# define neumann_pressure(i) (alpha.n[i] ? a.n[i]*fm.n[i]/alpha.n[i] :	\
			      a.n[i]*rho[]/(cm[] + SEPS))
#else
# define neumann_pressure(i) (a.n[i]*fm.n[i]/alpha.n[i])
#endif

p[right] = neumann (neumann_pressure(ghost));
p[left]  = neumann (- neumann_pressure(0));

#if AXI
uf.n[bottom] = 0.;
uf.t[bottom] = dirichlet(0); // since uf is multiplied by the metric which
                             // is zero on the axis of symmetry
p[top]    = neumann (neumann_pressure(ghost));
#else // !AXI
#  if dimension > 1
p[top]    = neumann (neumann_pressure(ghost));
p[bottom] = neumann (- neumann_pressure(0));
#  endif
#  if dimension > 2
p[front]  = neumann (neumann_pressure(ghost));
p[back]   = neumann (- neumann_pressure(0));
#  endif
#endif // !AXI

/**
For [embedded boundaries on trees](/src/embed-tree.h), we need to
define the pressure gradient for prolongation of pressure close to
embedded boundaries. */

#if TREE && EMBED
void pressure_embed_gradient (Point point, scalar p, coord * g)
{
  foreach_dimension()
    g->x = rho[]/(cm[] + SEPS)*(a.x[] + a.x[1])/2.;
}
#endif // TREE && EMBED

/**
## Initial conditions */

event defaults (i = 0)
{

  CFL = 0.8;

  /**
  The pressures are never dumped. */

  p.nodump = pf.nodump = true;
  
  /**
  The default density field is set to unity (times the metric). */

  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  }
  else if (!is_constant(alpha.x)) {
    face vector alphav = alpha;
    foreach_face()
      alphav.x[] = fm.x[];
    boundary ((scalar *){alpha});
  }

  /**
  On trees, refinement of the face-centered velocity field needs to
  preserve the divergence-free condition. */

#if TREE
  uf.x.refine = refine_face_solenoidal;

  /**
  When using [embedded boundaries](/src/embed.h), the restriction and
  prolongation operators need to take the boundary into account. */

#if EMBED
  uf.x.refine = refine_face;
  foreach_dimension()
    uf.x.prolongation = refine_embed_face_x;
  for (scalar s in {p, pf, u, g}) {
    s.restriction = restriction_embed_linear;
    s.refine = s.prolongation = refine_embed_linear;
  }
  for (scalar s in {p, pf})
    s.embed_gradient = pressure_embed_gradient;
#endif // EMBED
#endif // TREE
}

/**
After user initialisation, we initialise the face velocity and fluid
properties. */

/**
## Approximate projection

This function constructs the centered pressure gradient and
acceleration field *g* using the face-centered acceleration field *a*
and the cell-centered pressure field *p*. */

void centered_gradient (scalar deltap, vector deltag)
{

  /**
  We first compute a face field $\mathbf{g}_f$ combining both
  acceleration and pressure gradient. */

  face vector deltagf[];
  foreach_face()
  deltagf.x[] = fm.x[]*(a.x[] - aold.x[]) - alpha.x[]*(deltap[] - deltap[-1])/Delta;
//    gf.x[] = fm.x[]*a.x[] - alpha.x[]*(p[] - p[-1])/Delta;
  boundary_flux ({deltagf});

  /**
  We average these face values to obtain the centered, combined
  acceleration and pressure gradient field. */

  trash ({deltag});
  foreach()
  foreach_dimension()
  deltag.x[] = (deltagf.x[] + deltagf.x[1])/(fm.x[] + fm.x[1] + SEPS); //!!! ???: Weugene, numerator = gf.x[]*fm.x[] + gf.x[1]*fm.x[1]
  boundary ((scalar *){deltag});
}

double dtmax;

event init (i = 0)
{
  boundary ((scalar *){u});
  trash ({uf});
  foreach_face()
    uf.x[] = fm.x[]*face_value (u.x, 0);
  boundary ((scalar *){uf});
  mgp = project (uf, p, alpha, dt, mgp.nrelax);

  centered_gradient (p, g); //calc gf^{n+1}, g^{n+1} using p^{n+1}, a^{n+1/2}
  /**
  We update fluid properties. */

  event ("properties");

  /**
  We set the initial timestep (this is useful only when restoring from
  a previous run). */

  dtmax = DT;
  event ("stability");
}

/**
## Time integration

The timestep for this iteration is controlled by the CFL condition,
applied to the face centered velocity field $\mathbf{u}_f$; and the
timing of upcoming events. */

event set_dtmax (i++,last) dtmax = DT;

event stability (i++,last) {
  dt = dtnext (stokes ? dtmax : timestep (uf, dtmax));

  fprintf(ferr, "stability: dt=%g DT=%g dtmax=%g\n", dt, DT, dtmax );
}

/**
If we are using VOF or diffuse tracers, we need to advance them (to
time $t+\Delta t/2$) here. Note that this assumes that tracer fields
are defined at time $t-\Delta t/2$ i.e. are lagging the
velocity/pressure fields by half a timestep. */

event vof (i++,last);
event tracer_advection (i++,last);

/**
The fluid properties such as specific volume (fields $\alpha$ and
$\alpha_c$) or dynamic viscosity (face field $\mu_f$) -- at time
$t+\Delta t/2$ -- can be defined by overloading this event. */

event properties (i++,last) {
  boundary ({alpha, mu, rho, kappa}); // Weugene: kappa added
}

event tracer_diffusion (i++,last);

/**
### Predicted face velocity field

For second-order in time integration of the velocity advection term
$\nabla\cdot(\mathbf{u}\otimes\mathbf{u})$, we need to define the face
velocity field $\mathbf{u}_f$ at time $t+\Delta t/2$. We use a version
of the Bell-Collela-Glaz [advection scheme](/src/bcg.h) and the
pressure gradient and acceleration terms at time $t$ (stored in vector
$\mathbf{g}$). */

void prediction()
{
  vector du;
  foreach_dimension(){
    scalar s = new scalar;
    du.x = s;
  }

  if (u.x.gradient)
    foreach()
      foreach_dimension() {
#if EMBED
        if (!fs.x[] || !fs.x[1])
            du.x[] = 0.;
	    else
#endif
          du.x[] = u.x.gradient (u.x[-1], u.x[], u.x[1])/Delta;
      }
  else
    foreach()
        foreach_dimension() {
#if EMBED
            if (!fs.x[] || !fs.x[1])
	            du.x[] = 0.;
	        else
#endif
	            du.x[] = (u.x[1] - u.x[-1])/(2.*Delta);
    }
  boundary ((scalar *){du});

  trash ({uf});
  foreach_face() {
    double un = dt*(u.x[] + u.x[-1])/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    uf.x[] = u.x[i] + (g.x[] + g.x[-1])*dt/4. + s*(1. - s*un)*du.x[i]*Delta/2.;
    #if dimension > 1
    if (fm.y[i,0] && fm.y[i,1]) {
      double fyy = u.y[i] < 0. ? u.x[i,1] - u.x[i] : u.x[i] - u.x[i,-1];
      uf.x[] -= dt*u.y[i]*fyy/(2.*Delta);
    }
    #endif
    #if dimension > 2
    if (fm.z[i,0,0] && fm.z[i,0,1]) {
      double fzz = u.z[i] < 0. ? u.x[i,0,1] - u.x[i] : u.x[i] - u.x[i,0,-1];
      uf.x[] -= dt*u.z[i]*fzz/(2.*Delta);
    }
    #endif
    uf.x[] *= fm.x[];
  }
  boundary ((scalar *){uf});

  delete ((scalar *){du});
}

/**
### Advection term

We predict the face velocity field $\mathbf{u}_f$ at time $t+\Delta
t/2$ then project it to make it divergence-free. We can then use it to
compute the velocity advection term, using the standard
Bell-Collela-Glaz advection scheme for each component of the velocity
field. */

event advection_term (i++,last)
{
//  event("vtk_file");
  if (!stokes) {
    prediction();
#if BRINKMAN_PENALIZATION && MODIFIED_CHORIN
    mgpf = project_bp (uf, pf, alpha, 0.5*dt, mgpf.nrelax, fs, target_U, u, eta_chorin);
#else
    mgpf = project (uf, pf, alpha, dt/2., mgpf.nrelax);
#endif
    advection ((scalar *){u}, uf, dt, (scalar *){g});// original version
//    event("vtk_file");
  }
}

/**
### Viscous term

We first define a function which adds the pressure gradient and
acceleration terms. */

static void correction (double dt)
{
  foreach()
    foreach_dimension()
      u.x[] += dt*g.x[]; //original version
  boundary ((scalar *){u});
}

/**
The viscous term is computed implicitly. We first add the pressure
gradient and acceleration terms, as computed at time $t$, then call
the implicit viscosity solver. We then remove the acceleration and
pressure gradient terms as they will be replaced by their values at
time $t+\Delta t$. */

event viscous_term (i++,last)
{
  if (constant(mu.x) != 0.) {
    correction (dt);
    mgu = viscosity (u, mu, rho, dt, mgu.nrelax);
//    correction (-dt);//original
  }
//    event("vtk_file");
  /**
  We reset the acceleration field (if it is not a constant). */

  if (!is_constant(a.x)) {
    face vector af = a;
    trash ({af});
    foreach_face(){
      af.x[] = 0.;
    }
  }
}



/**
### Acceleration term

The acceleration term $\mathbf{a}$ needs careful treatment as many
equilibrium solutions depend on exact balance between the acceleration
term and the pressure gradient: for example Laplace's balance for
surface tension or hydrostatic pressure in the presence of gravity.

To ensure a consistent discretisation, the acceleration term is
defined on faces as are pressure gradients and the centered combined
acceleration and pressure gradient term $\mathbf{g}$ is obtained by
averaging. 

The (provisionary) face velocity field at time $t+\Delta t$ is
obtained by interpolation from the centered velocity field. The
acceleration term is added. */

event acceleration (i++,last)
{
  trash ({uf});
  boundary ((scalar *){a});
  foreach_face(){
    uf.x[] = fm.x[]*(face_value (u.x, 0) + dt*(a.x[] - aold.x[]));
  }
  boundary ((scalar *){uf});
//  event("vtk_file");
}


/**
To get the pressure field at time $t + \Delta t$ we project the face
velocity field (which will also be used for tracer advection at the
next timestep). Then compute the centered gradient field *g*. */

event projection (i++,last)
{
//  scalar deltap[];
//  vector deltag[];
//  trash({deltap, deltag});

//  foreach() {
//      ptotal[] = p[];
//      p[] = 0;
//  }
//      foreach(){
//        deltap[]=0;
//        foreach_dimension() deltag.x[]=0;
//      }
#if BRINKMAN_PENALIZATION && MODIFIED_CHORIN
  mgp = project_bp (uf, deltap, alpha, dt, mgp.nrelax, fs, target_U, u, eta_chorin);
#else
  mgp = project (uf, deltap, alpha, dt, mgp.nrelax);
#endif

//  event("vtk_file");
  centered_gradient (deltap, deltag); //calc gf^{n+1}, g^{n+1} using p^{n+1}, a^{n+1/2}
  /**
  We add the gradient field *g* to the centered velocity field. */
  foreach(){
    foreach_dimension()
    {
      u.x[] += dt*deltag.x[];
      g.x[] += deltag.x[];
    }
    p[] += deltap[];
  }
  boundary ((scalar *){u, g, p});

  if (!is_constant(a.x)) {
      face vector af = aold;
      foreach_face()
      {
          af.x[] = aold.x[];
      }
  }

//  event("vtk_file");
//  correction (dt);
}

//#if BRINKMAN_PENALIZATION
//event brinkman_penalization(i++, last){
//    brinkman_correction(u, uf, rho, dt);
//}
//#endif
/**
Some derived solvers need to hook themselves at the end of the
timestep. */

event end_timestep (i++, last);


/**
Output vtk files*/
event vtk_file (i++, last);// correct. Added by Weugene
/**
## Adaptivity

After mesh adaptation fluid properties need to be updated. When using
[embedded boundaries](/src/embed.h) the fluid fractions and face
fluxes need to be checked for inconsistencies. */

#if TREE
event adapt (i++,last) {
#if EMBED
  fractions_cleanup (cs, fs);
  foreach_face()
    // fs = 0 solid
    if (uf.x[] && !fs.x[])
      uf.x[] = 0.;
  boundary ((scalar *){uf});
#endif
//#if BRINKMAN_PENALIZATION
//  foreach_face()
//    // fs_face = 1 solid
//    if ((uf.x[] - target_Uf.x[]) && fs_face.x[])
//        uf.x[] = target_Uf.x[];
//  boundary ((scalar *){uf});
//#endif
  event ("properties");
}
#endif

/**
## See also

* [Double projection](double-projection.h)
* [Performance monitoring](perfs.h)
*/





#endif
