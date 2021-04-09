// Examples by Antoon
// https://groups.google.com/g/basilisk-fr/c/iJYpRYRHh9c
// A minimal .h module file contents:
#include "grid/octree.h"
#include "run.h"
// #define CONSTANT_V 1
(const) vector v;

//vector realv[];
event timestep (i++) {
  dt = dtnext(DT);
}

vector define_v(void); //prototype

event that_uses_v (i++) {
  v = define_v(); //May return const-typed or not
//  foreach(){
//    realv.x[] = noise();
//    realv.y[] = noise();
//    realv.z[] = noise();
//  }
//  foreach()
//    printf ("v: %g %g %g realv: %g %g %g\n", v.x[], v.y[], v.z[], realv.x[], realv.y[], realv.z[]);
  fprintf (ferr, "i=%d N=%d\n", i, N);
}

// A minimal .c setup file content:

// Two user defined v examples:
#if CONSTANT_V
  vector define_v (void) {
    double vx = t, vy = t*t, vz=2*t;
    const vector v[] = {vx, vy, vz};
    return v;
  }
#else
  vector define_v (void) {
    vector v[];
    foreach(){
      v.x[] = x*y*t*2;
      v.y[] = x*y*t*2;
      v.z[] = x*y*t*3;
      }
    return v;
  }
#endif

int main() {
  N = 1<<8;
  DT = 0.1;
  run();
}

event stop (i = 10000);