#include "grid/octree.h"
#include "run.h"
#include "timestep.h"
#include "utils.h"
#include "lambda2.h"
#include "../src_local/output_vtu_foreach.h"

#include <stdio.h>
#include <stdlib.h>
#include <wordexp.h>
#include <ctype.h>
scalar fs[], f[], p[];
vector u[];
double length_min = 1e+30, length_max = -1e+30, length = 1;
int i_take = 1;
double myt = 0, shiftm = 6, shiftp = 5;
wordexp_t fp;
char **w, dump_name[30];

//get time from dump name (dump-1.1234 => 1.1234)
double get_double(const char *str)
{
    /* First skip non-digit characters */
    /* Special case to handle negative numbers and the `+` sign */
    while (*str && !(isdigit(*str) || ((*str == '-' || *str == '+') && isdigit(*(str + 1)))))
        str++;
    /* The parse to a double */
    return strtod(str, NULL);
}

int main (int argc, char * argv[]) {
    // set which dump files will be converted: each $(i_take)th
    // by default each dump will be converted
    if (argc > 1)
        i_take = atoi (argv[1]);
    if (argc > 2)
        shiftm = fabs(atof (argv[2]));
    if (argc > 3)
        shiftp = fabs(atof (argv[3]));
    if (argc > 4)
        iter_fp = atoi (argv[4]);
    fprintf(ferr, "i_take=%d shiftm=%g shiftp=%g iter_fp=%d\n", i_take, shiftm, shiftp, iter_fp);
    wordexp("dump-*", &fp, 0);
    w = fp.we_wordv;
    for (int i = 0; i < fp.we_wordc; i++) fprintf(ferr, "All dump files: %s\n", w[i]);

    for (int i = 0; i < fp.we_wordc; i += i_take) {
        myt =  fabs(get_double(w[i]));
        strcpy(dump_name, w[i]);
        fprintf(ferr, "reading dump file: %s at time t= %g\n", dump_name, myt);
        run();
    }
    wordfree(&fp);
}

event init (t = 0) {
    bool success = restore (file = dump_name);
    fprintf(ferr, "file has been read: L0=%g\n", L0);

    if (!success) {
        fprintf(ferr, "can't open the file %s. Missing this file, go to the next file\n", dump_name);
        return 0;
    }
}

void calculate_aux_fields(vector u, scalar l, scalar omega, scalar l2){
    foreach() l[] = level;
    vorticity (u, omega);
    lambda2 (u, l2);
}

event vtk_file (i++)
{
    scalar l[], omega[], l2[];
    calculate_aux_fields(u, l, omega, l2);

    double xcg = 0, dvtmp, volume = 0, volumeg = 0 ;
    length_min = 1e+30, length_max = -1e+30, length = 0;
    foreach( reduction(+:xcg) reduction(+:volume) reduction(+:volumeg)) {
        if (fs[]<1){
        dvtmp = (1.0 - f[])*(1.0 - fs[])*dv(); // gas volume
        volumeg += dvtmp;//gas liquid
        volume += (1.0 - fs[])*dv();//channel volume
        xcg   += x*dvtmp;// Along x
        }
    }
    xcg /= volumeg;
    length_min = xcg - shiftm;
    length_max = xcg + shiftp;
    length = length_max - length_min;

    fprintf (ferr, "x= %g length_min= %g length_max= %g length= %g it_fp= %d\n"
                    "volume= %g volumeg= %g\n",
                    xcg, length_min, length_max, length, iter_fp,
                    volume, volumeg);
    if(0){
        char filename[80];
        sprintf(filename, "Conv_contour-%g.csv", t);
        FILE *fpc = fopen (filename, "w");
        fprintf(fpc, "x,y,z,f,u.x,u.y,u.z,omega,p\n");
        foreach(){
            if (f[]>0 && f[]<1) {
                fprintf(fpc, "%g,%g,%g,%g,%g,%g,%g,%g,%g\n",x,y,z,f[],u.x[],u.y[],u.z[],omega[],p[]);
            }
        }
        fclose (fpc);
    }else{
        char filename[80];
        sprintf(filename, "Conv_contour_short-%g.csv\n", t);
        FILE *fpc = fopen (filename, "w");
        fprintf(fpc, "x,y,z,f");
        foreach(){
            if (f[]>0 && f[]<1) {
                fprintf(fpc, "%g,%g,%g,%g\n",x,y,z,f[]);
            }
        }
        fclose (fpc);
    }
    unrefine ( (x < length_min || x > length_max) && level >= 1);
    unrefine ( (sq(y) + sq(z) > sq(0.55)) && level >= 1);
    
    char subname[80]; sprintf(subname, "dump2pvd_compressed");
    output_vtu_MPI( subname, myt, (scalar *) {fs, f, l, l2, omega, p}, (vector *) {u});
}

event stop(t = 100){ // t = 100 should  be sufficiently big in order to reach this event
    return 0;
};








//#include "grid/octree.h"
//#include "run.h"
//#include "timestep.h"
//#include "utils.h"
//#include "lambda2.h"
//#include "../src_local/output_vtu_foreach.h"
//
//#include <stdio.h>
//#include <stdlib.h>
//#include <wordexp.h>
//#include <ctype.h>
//scalar fs[], omega[], l2[], l[];
//scalar f[], * interfaces = {f};
//vector u[];
//double length_min = 1e+30, length_max = -1e+30, length = 1;
//int i_take = 1;
//double myt=0;
//wordexp_t fp;
//char **w, dump_name[30];
//
////get time from dump name (dump-1.1234 => 1.1234)
//double get_double(const char *str)
//{
//    /* First skip non-digit characters */
//    /* Special case to handle negative numbers and the `+` sign */
//    while (*str && !(isdigit(*str) || ((*str == '-' || *str == '+') && isdigit(*(str + 1)))))
//        str++;
//    /* The parse to a double */
//    return strtod(str, NULL);
//}
//
//int main (int argc, char * argv[]) {
//    // set which dump files will be converted: each $(i_take)th
//    // by default each dump will be converted
//    if (argc > 1)
//        i_take = atoi (argv[1]);
//
//    wordexp("dump-*", &fp, 0);
//    w = fp.we_wordv;
//    for (int i = 0; i < fp.we_wordc; i++) fprintf(ferr, "All dump files: %s\n", w[i]);
//
//    for (int i = 0; i < fp.we_wordc; i += i_take) {
//        myt =  fabs(get_double(w[i]));
//        strcpy(dump_name, w[i]);
//        fprintf(ferr, "reading dump file: %s at time t= %g\n", dump_name, myt);
//        run();
//    }
//    wordfree(&fp);
//}
//
//event init (t = 0) {
//    bool success = restore (file = dump_name);
//    fprintf(ferr, "file has been read: L0=%g\n", L0);
//
//    if (!success) {
//        fprintf(ferr, "can't open the file %s. Missing this file, go to the next file\n", dump_name);
//        return 0;
//    }
//}
//
//event calculate_aux_fields(i++){
//    foreach() l[] = level;
//    vorticity (u, omega);
//    lambda2 (u, l2);
//}
//event coarsen_grid(i++){
//    double xcg = 0, dvtmp, volume = 0, volumeg = 0 ;
//    length_min = 1e+30, length_max = -1e+30, length = 0;
//    foreach( reduction(+:xcg) reduction(+:volume) reduction(+:volumeg)) {
//        if (fs[]<1){
//            dvtmp = (1.0 - f[])*(1.0 - fs[])*dv(); // gas volume
//            volumeg += dvtmp;//gas liquid
//            volume += (1.0 - fs[])*dv();//channel volume
//            xcg   += x*dvtmp;// Along x
//        }
//    }
//    xcg /= volumeg;
//    length_min = xcg - 5;
//    length_max = xcg + 4;
//    length = length_max - length_min;
//
//    fprintf (ferr, "x= %g length_min= %g length_max= %g length= %g it_fp= %d\n"
//                   "volume= %g volumeg= %g\n",
//             xcg, length_min, length_max, length, iter_fp,
//             volume, volumeg);
//    unrefine ( (x < length_min || x > length_max) && level >= 1);
//    unrefine ( (sq(y) + sq(z) > sq(0.55)) && level >= 1);
//}
//
//
//
//event vtk_file (i++)
//{
//    char subname[80]; sprintf(subname, "dump2pvd_compressed");
//    output_vtu_MPI( subname, myt, (scalar *) {fs, f, l, l2, omega}, (vector *) {u});
//}
//
//event stop(t = 100){ // t = 100 should  be sufficiently big in order to reach this event
//    return 0;
//};


