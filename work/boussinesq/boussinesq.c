#include "convection_boussinesq.h"
#include "streamfunction.h"
//#include "navier-stokes/stream.h"
#include "output_fields/output_vtu_foreach.h"
#define MINLEVEL 4
#define MAXLEVEL 8

//scalar RT[];

double EndTime= 300.;



u.n[top] = dirichlet(0.);
u.n[bottom] = dirichlet(0.);
u.n[right] = dirichlet(0.);
u.n[left] = dirichlet(0.);

u.t[top] = dirichlet(0.);
u.t[bottom] = dirichlet(0.);
u.t[right] = dirichlet(0.);
u.t[left] = dirichlet(0.);



event init (t=0) {
    foreach(){
        T[] = x;
        foreach_dimension()
            u.x[] = 0.;
    }
    boundary ({T,u});
}


void backup_fields (scalar T, vector u, int nf)
{
    char name[80];
    FILE * fp ;
    nf > 0 ? sprintf(name, "cav2d_T_%6.6d.bin", nf) : sprintf(name, "cav2d_T.bin");
    fp = fopen(name, "w");
    output_matrix (T, fp, N, linear = true);
    fclose (fp);

    scalar psi[];
    streamfunction (u, psi);
    boundary ({psi});

    nf > 0 ? sprintf(name, "cav2d_Psi_%6.6d.bin", nf) : sprintf(name, "cav2d_Psi.bin");
    fp = fopen(name, "w");
    output_matrix (psi, fp, N, linear = true);
    fclose (fp);

    nf > 0 ? sprintf(name, "cav2d_%6.6d_n%3.3d.vtu", nf,pid()) : sprintf(name, "cav2d_n%3.3d.vtu",pid());
    fp = fopen(name, "w");
    output_vtu_ascii_foreach ((scalar *) {T,psi}, (vector *) {u}, N, fp, false);
    fclose (fp);
}

event logfile (t <= EndTime) {
    backup_fields(T,u,0);
}

scalar un[];
event init_un (i = 0) {
    foreach()
        un[] = u.x[];
}
#include "global_nusselt.h"
event logfile (t += 1.0; t <= EndTime) {
    static int nf = MINLEVEL;
    double deltau = change (u.x, un);
    if (deltau < 1e-8 && t > 5.0){
        backup_fields(T,u,N);

        FILE * fp ;
        if (nf > MINLEVEL) {
            fp = fopen("cav2d.asc", "a");
        } else {
            fp = fopen("cav2d.asc", "w");
            fprintf (fp, "[1]Ra [2]Pr [3]N [4]nuleft [5]nuright [6]nuvol [7]umin [8]umax [9]vmin [10]vmax\n");
        }

        double nu_l = 0., nu_r = 0, nu_vol=0. ;
        nu_r=nusselt_right(T); nu_l=nusselt_left(T); nu_vol=nusselt_vol(T,u);

        stats velx = statsf (u.x); stats vely = statsf (u.y);

        fprintf (fp, "%.9g %.9g %d %.9g %.9g %.9g %.9g %.9g %.9g %.9g \n",Ra,Pr,N,nu_l,nu_r,nu_vol, velx.min, velx.max, vely.min, vely.max );
        fclose (fp);

        refine (level < nf+1);
        if (nf >= MAXLEVEL) return 1;
        nf++; N*=2;
    }
  scalar div[];
  foreach() {
    div[] = 0.;
    foreach_dimension()
      div[] += u.x[1] - u.x[];
    div[] /= Delta;
  }
  stats s0 = statsf (div);
    fprintf (stderr, "%f %.9g %.9g %.9g \n", t, deltau, s0.sum/s0.volume, s0.max);
}

int main() {
    L0 = 1.;
    X0 = Y0 = -0.5;
    DT = 0.1;
    TOLERANCE = 1e-5;
    RT[left] = dirichlet(-0.5);
    T[right] = dirichlet(0.5);a = 1e3; Pr = 0.71; N = 1<<MINLEVEL ;
    run();
}




