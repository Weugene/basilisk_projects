#include "grid/quadtree.h"
#include "run.h"
#include "timestep.h"
#include "utils-weugene.h"
#include "output_htg.h"
scalar u[], uexact[];
double dtmax;
int LEVEL = 6;

int main(int argc, char * argv[]) {
    DT = 1e-4;
    CFL = 0.5;
    size(1.0);
    init_grid(1 << LEVEL);
    origin (-L0/2, -L0/2.);
    run();
    return 0;
}

#define u_BC (sin(3*pi*x)*sin(3*pi*y))
#define d2u_BC (-18*sq(pi)*u_BC)

//u[left] = neumann(0);
//u[right] = neumann(0);
//u[bottom] = neumann(0);
//u[top] = neumann (0);

u[left] = dirichlet(u_BC);
u[right] = dirichlet(u_BC);
u[bottom] = dirichlet(u_BC);
u[top] = dirichlet(u_BC);

//u[left] = (u_BC);
//u[right] = (u_BC);
//u[bottom] = (u_BC);
//u[top] = (u_BC);

//u[left] = 2*u_BC - u[];
//u[right] = 2*u_BC - u[];
//u[bottom] = 2*u_BC - u[];
//u[top] = 2*u_BC - u[];

event init (t = 0) {
    foreach(){
        u[] = 0;
    }
    event("report");
}

event set_dtmax (i++) dtmax = DT;

event stability (i++) {
    dt = dtnext ( dtmax );
}

event step(i++){
    u[left] = dirichlet(u_BC);
    u[right] = dirichlet(u_BC);
    u[bottom] = dirichlet(u_BC);
    u[top] = dirichlet(u_BC);
    foreach(){
        u[] = u[] + dt * ( (u[-1] - 4*u[] + u[1] + u[0,1] + u[0,-1])/sq(Delta) - d2u_BC);
    }
}

event report(t+=0.05){
    char path[] = "res"; // no slash at the end!!
    char prefix[] = "data";
    scalar l[];
    foreach() {
        l[] = level;
        uexact[] = u_BC;
    }
    output_htg(path, prefix,  (iter_fp) ? t + dt : 0, (scalar *) {u, uexact, l});
    if (i==10000 || i<=1){
        foreach() {
            double ubc = sin(3*pi*(x-Delta/2.))*sin(3*pi*(y));
            if (x < -L0 / 2. + L0 / pow(2., LEVEL))
                fprintf(ferr, "i=%d foreach left:xy= %g %g u[-1]= %12.10g u[]= %12.10g u[1]= %12.10g u_BC=%12.10g compare=%12.10g\n",
                    i, x, y, u[-1], u[], u[1], ubc, 2 * ubc - u[]);

        }
        foreach_boundary(left){
            fprintf(ferr, "i=%d foreach_boundary left:xy= %g %g u[-1]= %12.10g u[]= %12.10g u[1]= %12.10g u_BC=%12.10g compare=%12.10g\n",
                        i, x, y, u[-1], u[], u[1], u_BC, 2 * sin(3*pi*x)*sin(3*pi*y) - u[]);
        }
    }
}

event stop_end(t = 1){
}