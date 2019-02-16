/**
# Ohmic conduction flux of charged species

This function computes the fluxes due to ohmic conduction appearing in
the [Nernst--Planck
equation](http://en.wikipedia.org/wiki/Nernst%E2%80%93Planck_equation). The
species charge concentrations are then updated using the explicit
scheme
$$
c^{n+1}_i = c^n_i +\Delta t \, \nabla \cdot( K_i c^n_i \nabla \phi^n)
$$ 
where $c_i$ is the volume density of the $i$-specie, $K_i$ its volume
electric conductivity and $\phi$ the electric potential. */

extern scalar phi;

struct Species {
  scalar * c; // A list of the species concentration and their corresponding
  int * z;    // valences
  double dt;
  // optional
  vector * K; // electric mobility (default the valence)
};

void ohmic_flux (struct Species sp)
{
  /**
  If the volume conductivity is not provided it is set to the value of
  the valence. */
  
  if (!sp.K) { // fixme: this does not work yet
    int i = 0;
    for (scalar s in sp.c) {
      const face vector kc[] = {sp.z[i], sp.z[i]}; i++;
      sp.K = vectors_append (sp.K, kc);
    }
  }

  scalar c;
  (const) face vector K;
  for (c, K in sp.c, sp.K) {

    /**
    The fluxes of each specie through each face due to ohmic transport
    are */

    face vector f[];
    foreach_face()
      f.x[] = K.x[]*(c[] + c[-1,0])*(phi[] - phi[-1,0])/(2.*Delta);
    boundary_flux ({f});

    /**
    The specie concentration is updated using the net amount of that
    specie leaving/entering each cell through the face in the interval
    $dt$ */

    foreach()
      foreach_dimension()
        c[] += sp.dt*(f.x[1,0] - f.x[])/Delta;
    boundary ({c});    
  }
}
