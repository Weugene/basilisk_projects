/**
# An example of moving solid

This is a fairly rough way of approximating boundary conditions on a
moving (or stationary) solid.

A cylinder moves through a liquid with a velocity corresponding to a
Reynolds number of 160 i.e. this is equivalent in principle to the
[von Karman vortex street example](/src/examples/karman.c) but in a
different reference frame. */

#include "navier-stokes/centered.h"
#include "fractions.h"
#include "../src_local/output_vtu_foreach.h"
#include "tracer.h"
scalar f[];
scalar * tracers = {f};

f[right]    = dirichlet(y < 0);

int maxlevel = 12;
int minlevel = 4;
double rad = 0.0625;
scalar cylinder[], omega[], divu[];

int main(){
  L0 = 8.;
  DT=1e-2;
  origin (-0.5, -L0/2.);
  N = 512;
  const face vector muc[] = {0.00078125,0.00078125};
  mu = muc;
  run();
}

/**
We just make an (empty) long channel. */

event init (t = 0) {
    mask(y > 0.5 ? top: y < -0.5 ? bottom : none);
    if (!restore (file = "restart")) {
        int it = 0;
        do {
            fraction(f, (y <= 0) && (x > 0.5)? 1 : -1);
            fraction (cylinder, - (sq(x) + sq(y) - sq(rad)));
        }while (adapt_wavelet({f, cylinder}, (double []){1e-3, 1e-3}, maxlevel=maxlevel, minlevel=minlevel).nf != 0 && ++it <= 10);
        boundary(all);
    }

}

/**
The moving cylinder is created by first initialising a volume fraction
field at each timestep corresponding to the position and shape of the
moving cylinder. Note that *cylinder* should really be a temporary
field inside *moving_cylinder()*. We just keep it as a global field to
be able to display it in the movie below. */

event moving_cylinder (i++) {
  coord vc = {1.,0.}; // the velocity of the cylinder
  fraction (cylinder, - (sq(x - vc.x*t) + sq(y - vc.y*t) - sq(0.0625)));

  /**
  We then use this (solid) volume fraction field to impose the
  corresponding velocity field inside the solid (as a volume-weighted 
  average). */
  
  foreach()
    foreach_dimension()
      u.x[] = cylinder[]*vc.x + (1. - cylinder[])*u.x[];
  boundary ((scalar *){u});
}

/**
We output some convergence statistics... */

event logfile (i++){
    foreach() {
        divu[] = 0;
        foreach_dimension() divu[] += (uf.x[1]-uf.x[])/Delta;
    }
    double Linfu = -10;
    foreach( reduction(max:Linfu) ){
    if (fabs(divu[]) > Linfu) Linfu = fabs(divu[]);
    }
    fprintf (ferr, "i=%d t=%g dt=%g iter_p=%d iter_u=%d div u=%g \n", i, t, dt, mgp.i, mgu.i, Linfu);
}

/**
... and a movie of vorticity. 

![Animation of the vorticity field.](movingcylinder/vort.gif)
*/

event images (t += 0.1) {
  static FILE * fp = popen ("ppm2gif > vort.gif", "w");
  scalar omega[];
  vorticity (u, omega);

  /**
  Cells for which *m* is negative will be black in the movie. */

  scalar m[];
  foreach()
    m[] = 0.5 - cylinder[];
  boundary ({m});

  output_ppm (omega, fp, box = {{-0.5,-0.5},{7.5,0.5}}, mask = m,
	      min=-10, max=10, linear=true);
}

event vtk_file (t += 0.01){
    char subname[80]; sprintf(subname, "rk");
    scalar l[];
    scalar omega[], divu[];
    foreach() {
        divu[] = 0;
        foreach_dimension() divu[] += (uf.x[1]-uf.x[])/Delta;
    }
    vorticity (u, omega);
    foreach() {l[] = level; omega[] *= 1 - cylinder[];}
    output_vtu_MPI( (scalar *) {cylinder, f, omega, p, l, divu}, (vector *) {u}, subname, 0 );
}
/**
We adapt according to the error on the velocity field as usual. */

event adapt (i++) {
  adapt_wavelet ((scalar *){f, u}, (double[]){1e-3, 3e-2, 3e-2}, maxlevel=maxlevel, minlevel=minlevel);
}

event stop(t = 7);
/**
## See also

This method was used in Gerris in particular for these papers:

~~~bib
@InProceedings{wu2007,
  author =  {C. J. Wu and L. Wang},
  title =  {Direct numerical simulation of self-propelled 
swimming of 3D bionic fish school},
  booktitle =  {Computational Mechanics, Proceedings of ISCM 2007},
  year =  {2007},
  pages =  {1-17},
  url =  {http://gfs.sf.net/Wu-Wang-2007-Comp_Mech.pdf},
}

@Article{Lin-Lin2016,
author={Lin-Lin, Zhu and Hui, Guan and Chui-Jie, Wu},
title={Three-dimensional numerical simulation of a bird 
model in unsteady flight},
journal={Computational Mechanics},
year={2016},
pages={1--11},
issn={1432-0924},
doi={10.1007/s00466-015-1233-3},
url={http://dx.doi.org/10.1007/s00466-015-1233-3}
}
~~~
*/
