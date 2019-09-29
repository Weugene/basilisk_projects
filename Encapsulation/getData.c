/* Title: Janus Drops
# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
*/
#include "navier-stokes/centered.h"

scalar f1[], f2[];
char filename[80];
int nx,ny;
double xmin, ymin, xmax, ymax;
scalar * list = NULL;
scalar vel[], omega[];
event init(t = 0)
{
  restore (file = filename);
  boundary(all);
  foreach(){
    vel[] = sqrt(sq(u.x[]) + sq(u.y[]));
  }
  vorticity (u, omega);
  boundary (all);
  FILE * fp = ferr;
  nx++;
  ny++;
  double Deltax = 0.999999*(xmax-xmin)/(nx - 1);
  double Deltay = 0.999999*(ymax-ymin)/(ny - 1);
  int len = list_len(list);
  double ** field = (double **) matrix_new (nx, ny, len*sizeof(double));
  for (int i = 0; i < nx; i++) {
    double x = Deltax*i + xmin;
    for (int j = 0; j < ny; j++) {
      double y = Deltay*j + ymin;
      int k = 0;
      for (scalar s in list){
        field[i][len*j + k++] = interpolate (s, x, y);
      }
    }
  }

  for (int i = 0; i < nx; i++) {
    double x = Deltax*i + xmin;
    for (int j = 0; j < ny; j++) {
      double y = Deltay*j + ymin;
      fprintf (fp, "%g %g", x, y);
      int k = 0;
      for (scalar s in list){
        fprintf (fp, " %g", field[i][len*j + k++]);
      }
      fputc ('\n', fp);
    }
  }
  fflush (fp);
  fclose (fp);
  matrix_free (field);
}

int main(int a, char const *arguments[])
{
  sprintf (filename, "%s", arguments[1]);
  xmin = atof(arguments[2]);
  xmax = atof(arguments[3]);
  ymin = atof(arguments[4]);
  ymax = atof(arguments[5]);
  nx = atoi(arguments[6]);
  ny = atoi(arguments[7]);
  list = list_add (list, f1);
  list = list_add (list, f2);
  list = list_add (list, vel);
  list = list_add (list, omega);
  run();
}
