struct OutputVTK {
  scalar * list;
  FILE * fp;
  double shift;
};
#if dimension == 1
#define MY_BOX_CONDITION (x >= Pmin.x) && (x <= Pmax.x)
    #define MY_DELTA_BOX_CONDITION (x - 0.5*Delta >= Pmin.x) && (x + 0.5*Delta <= Pmax.x)
#elif dimension == 2
#define MY_BOX_CONDITION (x >= Pmin.x) && (x <= Pmax.x) && (y >= Pmin.y) && (y <= Pmax.y)
    #define MY_DELTA_BOX_CONDITION (x - 0.5*Delta >= Pmin.x) && (x + 0.5*Delta <= Pmax.x) && (y - 0.5*Delta >= Pmin.y) && (y + 0.5*Delta <= Pmax.y)
#elif dimension > 2
#define MY_BOX_CONDITION (x >= Pmin.x) && (x <= Pmax.x) && (y >= Pmin.y) && (y <= Pmax.y) && (z >= Pmin.z) && (z <= Pmax.z)
    #define MY_DELTA_BOX_CONDITION (x - 0.5*Delta >= Pmin.x) && (x + 0.5*Delta <= Pmax.x) && (y - 0.5*Delta >= Pmin.y) && (y + 0.5*Delta <= Pmax.y) && (z - 0.5*Delta >= Pmin.z) && (z + 0.5*Delta <= Pmax.z)
#endif
void output_vtk (struct OutputVTK ov)
{
  scalar * list = ov.list;
  FILE * fp = ov.fp;
  double shift = ov.shift;
  coord Pmin = {X0 + shift, Y0 + shift, Z0 + shift};
  coord Pmax = {X0 + L0 - shift, Y0 + L0 - shift, Z0 + L0 - shift};

  fputs ("# vtk DataFile Version 2.0\n"
	 "Basilisk\n"
	 "ASCII\n"
	 "DATASET UNSTRUCTURED_GRID\n \n", fp);

  vertex scalar psi[];
  int np = 0;
  foreach_vertex () {
     if (MY_BOX_CONDITION) {
     	psi[] = np;
     	np++;
     }
//     else{fprintf(ferr, "x=%g y=%g\n", x, y);}
  }
//    fprintf(ferr, "np = %d\n", np);
//  int npt=0;
//  foreach_vertex () {
//        npt++;
//  }
//  fprintf(ferr, "npt = %d\n", npt);
  fprintf (fp, "POINTS %d double\n", np);

  foreach_vertex () {
     if (MY_BOX_CONDITION)
#if dimension == 1
       fprintf (fp, "%g 0 0\n", x, y);
#elif dimension == 2
	   fprintf (fp, "%g %g 0\n", x, y);
#elif dimension > 2
       fprintf (fp, "%g %g %g\n", x, y, z);
#endif
  }

  fprintf (fp, "\n");

  int ncells=0;
  foreach () {
     if (MY_DELTA_BOX_CONDITION) ncells++;
  }

  fprintf (fp, "CELLS %i %i \n", ncells, ncells*5);

  foreach () {
    if (MY_DELTA_BOX_CONDITION)
#if dimension == 1
    fprintf (fp, "2 %i %i %i %i \n", (int) psi[], (int) psi[1]);
#elif dimension == 2
     fprintf (fp, "4 %i %i %i %i \n", (int) psi[], (int) psi[0,1], (int) psi[1,0], (int) psi[1,1]);
#else dimension > 2
     fprintf (fp, "8 %i %i %i %i \n", (int) psi[], (int) psi[0,1,0], (int) psi[1,0,0], (int) psi[1,1,0], (int) psi[0,0,1], (int) psi[0,1,1], (int) psi[1,0,1], (int) psi[1,1,1]);
#endif
  }
  fprintf (fp, "\n");


  fprintf (fp, "CELL_TYPES %i \n", ncells);
  for (int i = 0; i < ncells; i++)
	fprintf (fp, "8 \n");
  fprintf (fp, "\n");

  fprintf (fp, "POINT_DATA %d\n", np);
  for (scalar s in list) {
    fprintf (fp, "SCALARS %s double\n", s.name);
    fputs ("LOOKUP_TABLE default\n", fp);
    foreach_vertex () {
        if (MY_BOX_CONDITION) fprintf(fp, "%g\n", s[]);
      }
    fprintf (fp, "\n");
  }

  fflush (fp);
}