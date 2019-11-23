#ifndef BASILISK_HEADER_8
#define BASILISK_HEADER_8
#line 1 "./../src_local/penalization.h"
struct Penalization {
    vector u;
    double dt;
};
#ifdef BRINKMAN_PENALIZATION

    extern scalar fs;
    double eta_s = 1e-3, nu_s = 0, lambda_slip = 0;
    (const) vector n_sol = zerof, target_U = zerof;
    (const) face vector n_sol_face = zerof, target_U_face = zerof;
//    (const) scalar a_br = unity;
//    (const) scalar b_br = unity;
#ifdef DEBUG_BRINKMAN_PENALIZATION
    vector dbp[];
#endif
#endif

event brinkman_penalization(i++, last){
  vector utau[];
  double ubyn, d, velo;
  foreach() {
    ubyn = 0;
    foreach_dimension() ubyn += (u.x[] - target_U.x[])*n_sol.x[];
    foreach_dimension() utau.x[] = u.x[] - ubyn*n_sol.x[];
  }
  foreach_dimension(){
    vector grad_u[];
    gradients({utau.x}, {grad_u});
    foreach() {
      d = 0;
      foreach_dimension() d += grad_u.x[]*n_sol.x[];
      velo = target_U.x[] + lambda_slip*d;
      u.x[] = (u.x[] + (fs[]*dt/eta_s)*velo)/(1 + fs[]*dt/eta_s);
#ifdef DEBUG_BRINKMAN_PENALIZATION
      dbp.x[] = - (fs[]/eta_s)*(u.x[] - velo);
#endif
    }
  }

  double fs_face;
  double uftau, grad_uf;
  foreach_face(){
    fs_face = 0.5*(fs[-1] + fs[]);
    d = 0; uftau = 0;
    foreach_dimension() d += uf.x[]*n_sol_face.x[];
    foreach_dimension() uftau = uf.x[] - d*n_sol_face.x[];
    grad_uf = (uf.x[1] - uf.x[])*n_sol_face.x[]/Delta
#if dimension > 1
    + (uf.x[0,1] - uf.x[])*n_sol_face.y[]/Delta
#endif
#if dimension > 2
    + (uf.x[0,0,1] - uf.x[])*n_sol_face.z[]/Delta
#endif
    ;
    uf.x[] = (uf.x[] + (fs_face*dt/eta_s)*(target_U_face.x[] + lambda_slip*grad_uf))/(1 +fs_face*dt/eta_s)
  }
}
#endif
