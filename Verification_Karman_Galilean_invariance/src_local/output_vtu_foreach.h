#ifndef BASILISK_HEADER_3
#define BASILISK_HEADER_3
#line 1 "./../src_local/output_vtu_foreach.h"
/*
This function writes one XML file which allows to read the *.vtu files generated
by output_vtu_ascii_foreach() when used in MPI. Tested in (quad- and oct-)trees
using MPI.
*/
#define SMALL_VAL 1e-12

#ifdef _MPI
#if dimension == 1
    #define MY_BOX_CONDITION (!periodic_bc) || (x - Pmin.x - 0.5*Delta > 0) && (Pmax.x - x - 0.5*Delta > 0)
#elif dimension == 2
	#define MY_BOX_CONDITION (!periodic_bc) || (x - Pmin.x - 0.5*Delta > 0) && (Pmax.x - x - 0.5*Delta > 0) && (y - Pmin.y - 0.5*Delta > 0) && (Pmax.y - y - 0.5*Delta > 0)
#elif dimension > 2
    #define MY_BOX_CONDITION (!periodic_bc) || (x - Pmin.x - 0.5*Delta > 0) && (Pmax.x - x - 0.5*Delta > 0) && (y - Pmin.y - 0.5*Delta > 0) && (Pmax.y - y - 0.5*Delta > 0) && (z - Pmin.z - 0.5*Delta > 0) && (Pmax.z - z - 0.5*Delta > 0)
#endif
#else
#if dimension == 1
    #define MY_BOX_CONDITION (!periodic_bc) || (Pmax.x - x - 0.5*Delta > 0)
#elif dimension == 2
    #define MY_BOX_CONDITION (!periodic_bc) || (Pmax.x - x - 0.5*Delta > 0) && (Pmax.y - y - 0.5*Delta > 0)
#elif dimension > 2
    #define MY_BOX_CONDITION (!periodic_bc) || (Pmax.x - x - 0.5*Delta > 0) && (Pmax.y - y - 0.5*Delta > 0) && (Pmax.z - z - 0.5*Delta > 0)
#endif
#endif

@include <sys/types.h>
@include <sys/stat.h>
@include <unistd.h>
void output_pvtu_ascii (scalar * list, vector * vlist, int n, FILE * fp, char * subname)
{
    fputs ("<?xml version=\"1.0\"?>\n"
    "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
    fputs ("\t <PUnstructuredGrid GhostLevel=\"0\">\n", fp);
    fputs ("\t\t\t <PCellData Scalars=\"scalars\">\n", fp);
    for (scalar s in list) {
      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
      fputs ("\t\t\t\t </PDataArray>\n", fp);
    }
    for (vector v in vlist) {
      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n", v.x.name);
      fputs ("\t\t\t\t </PDataArray>\n", fp);
    }
    fputs ("\t\t\t </PCellData>\n", fp);
    fputs ("\t\t\t <PPoints>\n", fp);
    fputs ("\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
    fputs ("\t\t\t\t </PDataArray>\n", fp);
    fputs ("\t\t\t </PPoints>\n", fp);

    for (int i = 0; i < npe(); i++)
      fprintf (fp, "<Piece Source=\"%s_n%3.3d.vtu\"/> \n", subname, i);

    fputs ("\t </PUnstructuredGrid>\n", fp);
    fputs ("</VTKFile>\n", fp);
}

/*
This function writes one XML VTK file per PID process of type unstructured grid
(*.vtu) which can be read using Paraview. File stores scalar and vector fields
defined at the center points. Results are recorded on ASCII format. If one writes
one *.vtu file per PID process this function may be combined with
output_pvtu_ascii() above to read in parallel. Tested in (quad- and oct-)trees
using MPI. Also works with solids (when not using MPI).
*/
void output_vtu_ascii_foreach (scalar * list, vector * vlist, int n, FILE * fp, bool linear, double shift)
{
    int dim = 3; bool periodic_bc = shift > 0;
    coord Pmin = {X0 + SMALL_VAL, Y0 + SMALL_VAL, Z0 + SMALL_VAL};
    coord Pmax = {X0 + L0 - SMALL_VAL, Y0 + L0 - SMALL_VAL, Z0 + L0 - SMALL_VAL};
#if defined(_OPENMP)
  int num_omp = omp_get_max_threads();
  omp_set_num_threads(1);
#endif

  vertex scalar marker[];
  int no_points = 0, no_cells=0 ;
  foreach_vertex(){
    if (MY_BOX_CONDITION) {
      marker[] = no_points++;
    }else{
	  marker[] = -1;
    }
  }
  foreach(){
    if (MY_BOX_CONDITION) no_cells += 1;
  }

  fputs ("<?xml version=\"1.0\"?>\n"
  "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fputs ("\t <UnstructuredGrid>\n", fp);
  fprintf (fp,"\t\t <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", no_points, no_cells);
  fputs ("\t\t\t <CellData Scalars=\"scalars\">\n", fp);
  for (scalar s in list) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"ascii\">\n", s.name);
    foreach(){
      if (MY_BOX_CONDITION)
        fprintf (fp, "\t\t\t\t\t %g\n", val(s));
    }
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  for (vector v in vlist) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" Name=\"%s\" format=\"ascii\">\n", v.x.name);
    foreach(){
      if (MY_BOX_CONDITION)
      #if dimension == 1
        fprintf (fp, "\t\t\t\t\t %g 0 0.\n", val(v.x));
      #endif
      #if dimension == 2
          fprintf (fp, "\t\t\t\t\t %g %g 0.\n", val(v.x), val(v.y));
      #endif
      #if dimension > 2
          fprintf (fp, "\t\t\t\t\t %g %g %g\n", val(v.x), val(v.y), val(v.z));
      #endif
    }
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  fputs ("\t\t\t </CellData>\n", fp);
  fputs ("\t\t\t <Points>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n", fp);
  foreach_vertex(){
    if (MY_BOX_CONDITION)
    #if dimension == 1
      fprintf (fp, "\t\t\t\t\t %g %g %g\n", x, y, z);
    #endif
    #if dimension == 2
      fprintf (fp, "\t\t\t\t\t %g %g %g\n", x, y, z);
    #endif
    #if dimension > 2
      fprintf (fp, "\t\t\t\t\t %g %g %g\n", x, y, z);
    #endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Points>\n", fp);
  fputs ("\t\t\t <Cells>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n", fp);
  foreach(){
    if (MY_BOX_CONDITION)
    #if dimension == 1
      fprintf (fp, "\t\t\t\t\t %d %d \n", (int)marker[], (int)marker[1]);
    #endif
    #if dimension == 2
      fprintf (fp, "\t\t\t\t\t %d %d %d %d \n", (int)marker[], (int)marker[1,0], (int)marker[1,1], (int)marker[0,1]);
    #endif
    #if dimension > 2
      fprintf (fp, "\t\t\t\t\t %d %d %d %d %d %d %d %d \n", (int)marker[], (int)marker[1,0,0], (int)marker[1,1,0], (int)marker[0,1,0], (int)marker[0,0,1], (int)marker[1,0,1], (int)marker[1,1,1], (int)marker[0,1,1]);
    #endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n", fp);

  for (int i = 1; i < no_cells+1; i++){
#if dimension == 2
    fprintf (fp, "\t\t\t\t\t %d \n", i*4);
#endif
#if dimension > 2
    fprintf (fp, "\t\t\t\t\t %d \n", i*8);
#endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", fp);
  foreach(){
    if (MY_BOX_CONDITION)
    #if dimension == 1
      fputs ("\t\t\t\t\t 3 \n", fp); //VTK_LINE (=3)
    #endif
    #if dimension == 2
      fputs ("\t\t\t\t\t 9 \n", fp); //VTK_QUAD (=9)
    #endif
    #if dimension > 2
      fputs ("\t\t\t\t\t 12 \n", fp); //VTK_HEXAHEDRON
    #endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Cells>\n", fp);
  fputs ("\t\t </Piece>\n", fp);
  fputs ("\t </UnstructuredGrid>\n", fp);
  fputs ("</VTKFile>\n", fp);
  fflush (fp);
#if defined(_OPENMP)
  omp_set_num_threads(num_omp);
#endif
}

/*
This function writes one XML file which allows to read the *.vtu files generated
by output_vtu_bin_foreach() when used in MPI. Tested in (quad- and oct-)trees
using MPI.
*/
void output_pvtu_bin (scalar * list, vector * vlist, int n, FILE * fp, char * subname)
{
    int dim = 3;
    fputs ("<?xml version=\"1.0\"?>\n"
    "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
    fputs ("\t <PUnstructuredGrid GhostLevel=\"0\">\n", fp);
    fputs ("\t\t\t <PCellData Scalars=\"scalars\">\n", fp);
    for (scalar s in list) {
      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" Name=\"%s\" format=\"appended\"/>\n", s.name);
    }
    for (vector v in vlist) {
      fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"%d\" Name=\"%s\" format=\"appended\"/>\n", dim, v.x.name);
    }
    fputs ("\t\t\t </PCellData>\n", fp);
    fputs ("\t\t\t <PPoints>\n", fp);
    fprintf (fp,"\t\t\t\t <PDataArray type=\"Float64\" NumberOfComponents=\"%d\" format=\"ascii\"/>\n", dim);
    fputs ("\t\t\t </PPoints>\n", fp);

    for (int i = 0; i < npe(); i++)
      fprintf (fp, "<Piece Source=\"%s_n%3.3d.vtu\"/> \n", subname, i);

    fputs ("\t </PUnstructuredGrid>\n", fp);
    fputs ("</VTKFile>\n", fp);
}

/*
This function writes one XML VTK file per PID process of type unstructured grid
(*.vtu) which can be read using Paraview. File stores scalar and vector fields
defined at the center points. Results are recorded on binary format. If one writes
one *.vtu file per PID process this function may be combined with
output_pvtu_bin() above to read in parallel. Tested in (quad- and oct-)trees
using MPI. Also works with solids (when not using MPI).
*/

void output_vtu_bin_foreach (scalar * list, vector * vlist, int n, FILE * fp, bool linear, double shift)
{
  int dim = 3; bool periodic_bc = shift > 0;
  coord Pmin = {X0 + SMALL_VAL, Y0 + SMALL_VAL, Z0 + SMALL_VAL};
  coord Pmax = {X0 + L0 - SMALL_VAL, Y0 + L0 - SMALL_VAL, Z0 + L0 - SMALL_VAL};
#if defined(_OPENMP)
  int num_omp = omp_get_max_threads();
  omp_set_num_threads(1);
#endif
//  fprintf(fp1, "X0=%.12g Y0=%.12g L0=%.12g shift=%.12g Pmin.x=%.12g Pmin.y=%.12g Pmax.x=%.12g Pmax.y=%.12g\n", X0, Y0, L0, shift, Pmin.x, Pmin.y, Pmax.x, Pmax.y);

  vertex scalar marker[];
  int no_points = 0, no_cells = 0;
  foreach_vertex(){
    if (MY_BOX_CONDITION) {
      marker[] = no_points++;
    }else{
    	marker[] = -1; //if you see -1 in vtu file=> there is a mistake
    }
  }
  foreach(){
    if (MY_BOX_CONDITION) no_cells++;
  }
  fputs ("<?xml version=\"1.0\"?>\n"
  "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fputs ("\t <UnstructuredGrid>\n", fp);
  fprintf (fp,"\t\t <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", no_points, no_cells);
  fputs ("\t\t\t <CellData Scalars=\"scalars\">\n", fp);
  int count = 0;
  for (scalar s in list) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"appended\" offset=\"%d\">\n", s.name,count);
    count += ((no_cells)+1)*8;
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  for (vector v in vlist) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" NumberOfComponents=\"%d\"  format=\"appended\" offset=\"%d\">\n", v.x.name, dim, count);
    count += (no_cells*dim+1)*8;
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  fputs ("\t\t\t </CellData>\n", fp);
  fputs ("\t\t\t <Points>\n", fp);
  fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"%d\"  format=\"appended\" offset=\"%d\">\n", dim, count);
  count += (no_points*dim+1)*8;
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Points>\n", fp);
  fputs ("\t\t\t <Cells>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n", fp);
  foreach(){
    if (MY_BOX_CONDITION) {
#if dimension == 1
	    fprintf (fp, "\t\t\t\t\t %d %d \n", (int)marker[], (int)marker[1]);
#endif
#if dimension == 2
	    fprintf (fp, "\t\t\t\t\t %d %d %d %d \n", (int)marker[], (int)marker[1,0], (int)marker[1,1], (int)marker[0,1]);
#endif
#if dimension > 2
	    fprintf (fp, "\t\t\t\t\t %d %d %d %d %d %d %d %d \n", (int)marker[], (int)marker[1,0,0], (int)marker[1,1,0], (int)marker[0,1,0], (int)marker[0,0,1], (int)marker[1,0,1], (int)marker[1,1,1], (int)marker[0,1,1]);
#endif
    }
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n", fp);
  for (int i = 1; i < no_cells+1; i++){
#if dimension == 1
    fprintf (fp, "%d \n", i*2);
#endif
#if dimension == 2
    fprintf (fp, "%d \n", i*4);
#endif
#if dimension > 2
    fprintf (fp, "%d \n", i*8);
#endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", fp);
  foreach(){
    if (MY_BOX_CONDITION)
    #if dimension == 1
      fputs ("3 \n", fp); //VTK_LINE (=3)
    #endif
    #if dimension == 2
      fputs ("9 \n", fp); //VTK_QUAD (=9)
    #endif
    #if dimension > 2
      fputs ("12 \n", fp); //VTK_HEXAHEDRON
    #endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Cells>\n", fp);
  fputs ("\t\t </Piece>\n", fp);
  fputs ("\t </UnstructuredGrid>\n", fp);
  fputs ("\t <AppendedData encoding=\"raw\">\n", fp);
  fputs ("_", fp);
  unsigned long long block_len=no_cells*8;
#if dimension == 1
  double y=0, vy=0;
  double z=0, vz=0;
#endif
#if dimension == 2
  double z=0, vz=0;
#endif
  for (scalar s in list) {
    fwrite (&block_len, sizeof (unsigned long long), 1, fp);
    foreach()
      if (MY_BOX_CONDITION)
        fwrite (&val(s), sizeof (double), 1, fp);
  }
  block_len=no_cells*8*dim;
  for (vector v in vlist) {
    fwrite (&block_len, sizeof (unsigned long long), 1, fp);
    foreach(){
      if (MY_BOX_CONDITION){
      #if dimension == 1
        fwrite (&val(v.x), sizeof (double), 1, fp);
        fwrite (&vy, sizeof (double), 1, fp);
        fwrite (&vz, sizeof (double), 1, fp);
      #endif
      #if dimension == 2
        fwrite (&val(v.x), sizeof (double), 1, fp);
        fwrite (&val(v.y), sizeof (double), 1, fp);
        fwrite (&vz, sizeof (double), 1, fp);
      #endif
      #if dimension > 2
        fwrite (&val(v.x), sizeof (double), 1, fp);
        fwrite (&val(v.y), sizeof (double), 1, fp);
        fwrite (&val(v.z), sizeof (double), 1, fp);
      #endif
      }
    }
  }
  block_len=no_points*8*dim;
  fwrite (&block_len, sizeof (unsigned long long), 1, fp);
  foreach_vertex(){
    if (MY_BOX_CONDITION){
      fwrite (&x, sizeof (double), 1, fp);
      fwrite (&y, sizeof (double), 1, fp);
      fwrite (&z, sizeof (double), 1, fp);
    }
  }
  fputs ("\t\n", fp);
  fputs ("\t </AppendedData>\n", fp);
  fputs ("</VTKFile>\n", fp);
  fflush (fp);
#if defined(_OPENMP)
  omp_set_num_threads(num_omp);
#endif
}

void output_pvd_file(FILE * fp, int nf, float * file_timesteps, char * subname){
//    fputs("<?xml version=\"1.0\"?>\n",fp);
    fputs ("<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
           "\t <Collection>\n", fp);
    for (int i=0; i<=nf; i++) {
        fprintf(fp, "\t\t<DataSet timestep=\"%g\" part=\"0\" file=\"res/%s_0_%4.4d.pvtu\"/>\n", file_timesteps[i], subname, i);
    }//12.9g
    fputs ("\t </Collection>\n "
           "</VTKFile>\n", fp);
}
/* output_vtu_MPI produces *.pvtu files and *.vtu files. The user needs to specify list of scalars and vectors and subname.
*/
struct stat st = {0};
static int iter_fp=0;
static float file_timesteps[9999];
void output_vtu_MPI(scalar * list, vector * vlist, char * subname, double shift){
    int nf = iter_fp;
    char name_vtu[80];
    if (iter_fp == 0) {
        if (stat("res", &st) == -1) {
            mkdir("res", 0755);
        }
    }
    FILE *fp;
    if (nf>9999) { fprintf(stderr, "too many files, more than 9999"); exit(1); }
    sprintf(name_vtu, "res/%s_%4.4d_n%3.3d.vtu", subname, nf, pid());
    fp = fopen(name_vtu, "w");
    output_vtu_bin_foreach(list, vlist, 64, fp, true, shift);//64 and true is useless. It needs to support the interface
    fclose(fp);
    if (pid() == 0) {
        //pvtu file
        char name_pvtu[80], tmp[80];
	    sprintf(name_pvtu, "res/%s_0_%4.4d.pvtu", subname, nf);
        sprintf(tmp, "%s_%4.4d", subname, nf);
        fp = fopen(name_pvtu, "w");
        output_pvtu_bin(list, vlist, 64, fp, tmp);
        fclose(fp);
        //pvd file with timesteps
        char name_pvd[80];
        sprintf(name_pvtu, "%s.pvd", subname);
        fp = fopen(name_pvtu, "w");
        file_timesteps[nf] = t;
        output_pvd_file(fp, nf, file_timesteps, subname);
        fclose(fp);
    }
    @if _MPI
        MPI_Barrier(MPI_COMM_WORLD);
    @endif
#ifdef DEBUG_OUTPUT_VTU_MPI
    fprintf (ferr, "iter_fp: %d t=%g dt=%g\n", nf, t, dt);
#endif
    iter_fp++;
}


void face_vector2vector(face vector fv, vector mapped_data_lower, vector mapped_data_upper){
//    face vector kappa;
//    kappa = some_face_data;
//    vector mapped_data_lower, mapped_data_upper;
    foreach()
    foreach_dimension(){
        mapped_data_lower.x[] = fv.x[];
        mapped_data_upper.x[] = fv.x[1];
    }
}

#endif
