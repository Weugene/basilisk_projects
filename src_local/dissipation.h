//
// Created by Weugene on 01.04.2022.
//

struct Dissipation {
    scalar dis;
    vector u;
    face vector mu;
};

void dissipation (struct Dissipation p)
{
    vector u = p.u;
    scalar dis = p.dis;
    (const) face vector mu = p.mu;

#if TREE
    /* conservative coarse/fine discretisation (2nd order) */
    foreach () {
        dis[] = 0.;
    }

    foreach_dimension() {
        face vector taux[];
        foreach_face(x) {
            taux.x[] =  2.*mu.x[]*sq((u.x[] - u.x[-1])/Delta);
        }
        #if dimension > 1
        foreach_face(y)
            taux.y[] = mu.y[]*(u.x[] - u.x[0,-1] +
                   (u.y[1,-1] + u.y[1,0])/4. -
                   (u.y[-1,-1] + u.y[-1,0])/4.)
                  *(u.x[] - u.x[0,-1])/sq(Delta);
        #endif
        #if dimension > 2
        foreach_face(z)
            taux.z[] = mu.z[]*(u.x[] - u.x[0,0,-1] +
                   (u.z[1,0,-1] + u.z[1,0,0])/4. -
                   (u.z[-1,0,-1] + u.z[-1,0,0])/4.)
                  *(u.x[] - u.x[0,0,-1])/sq(Delta);
        #endif
        foreach () {
            double d = 0;
            foreach_dimension()
                d += (taux.x[1] + taux.x[]);
            dis[] += (d)/2.;
        }
    }
#else
/*   /\* "naive" discretisation (only 1st order on trees) *\/ */
    foreach () {
        dis[] = 0.;
        foreach_dimension()
            dis[] += ((mu.x[1,0]*sq(u.x[1] - u.x[])
                  + mu.x[]*sq(u.x[] - u.x[-1])
#if dimension > 1
                  + mu.y[0,1]*(u.x[0,1] - u.x[] +
		(u.y[1,0] + u.y[1,1])/4. -
		(u.y[-1,0] + u.y[-1,1])/4.)*(u.x[0,1] - u.x[])/2.
		+ mu.y[]*(u.x[] - u.x[0,-1] +
		(u.y[1,-1] + u.y[1,0])/4. -
		(u.y[-1,-1] + u.y[-1,0])/4.)*(u.x[] - u.x[0,-1])/2.
#endif
#if dimension > 2
                   + mu.z[0,0,1]*(u.x[0,0,1] - u.x[] +
        (u.z[1,0,0] + u.z[1,0,1])/4. -
        (u.z[-1,0,0] + u.z[-1,0,1])/4.)*(u.x[0,0,1] - u.x[])/2.
		- mu.z[]*(u.x[] - u.x[0,0,-1] +
        (u.z[1,0,-1] + u.z[1,0,0])/4. -
        (u.z[-1,0,-1] + u.z[-1,0,0])/4.)*(u.x[] - u.x[0,0,-1])/2.
#endif
      ) )/sq(Delta);
    }
#endif
}