//#include "run.h"
//
//// The "Solver" code
//scalar s[];
//
//event event1 (i++) {
//boundary({s});
//foreach_boundary(top){
//    fprintf(ferr, "event1: %g %g %g 2 \n", s[-1], s[], s[1]);
//    assert(s[ghost] == 2);
//}
//}
//
//event event2 (i++) {
//boundary({s});
//foreach_boundary(top){
//    fprintf(ferr, "event2: %g %g %g 1 \n", s[-1], s[], s[1]);
//    assert (s[ghost] == 1);
//}
//return 1;
//}
//// The setup code:
//int main() {
//    run();
//}
//// Hooks for the solver events..
//event event1 (i++) {
//    s[top] = 2;
//}
//
//event event2 (i++) {
//    s[top] = 1;
//}




#include "run.h"

// The "Solver" code
face vector uf[];
scalar s[];
event event1 (i++) {
  boundary({uf,s});

  foreach_boundary(left){
//          fprintf(ferr,"x=%g, y=%g: %g %g %g gh=%g %d\n", x, y, uf.x[-1], uf.x[], uf.x[1], uf.x[ghost], ghost );
          assert(uf.x[] == 22);
  }
  fprintf(ferr,"left\n" );

  foreach_boundary(right){
//          fprintf(ferr,"x=%g, y=%g: %g %g %g gh=%g %d\n", x, y, uf.x[-1], uf.x[], uf.x[1], uf.x[ghost], ghost );
          assert(uf.x[ghost] == 44);
  }
  fprintf(ferr,"right\n" );

  foreach_boundary(bottom){
//        fprintf(ferr,"x=%g, y=%g: %g %g %g\n", x, y, uf.y[-1], uf.y[], uf.y[1] );
        assert(uf.y[] == 2);
  }
  fprintf(ferr,"bottom\n" );

  foreach_boundary(top){
//      fprintf(ferr,"x=%g, y=%g: %g %g %g\n", x, y, uf.y[-1], uf.y[], uf.y[1] );
      assert(uf.y[ghost] == 4);
  }
  fprintf(ferr,"top\n" );

foreach_boundary(top){
//    fprintf(ferr, "event1: x=%g y=%g | %g %g %g %g \n", x, y, s[-1], s[], s[1], s[ghost]);
    assert(s[ghost] == y);
}

}

event event2 (i++) {
  boundary({s,uf});
  foreach_boundary(left){
          fprintf(ferr,"x=%g, y=%g: %g %g %g gh=%g %d\n", x, y, uf.x[-1], uf.x[], uf.x[1], uf.x[ghost], ghost );
          assert(uf.x[] == x);
  }
  fprintf(ferr,"left\n" );

  foreach_boundary(right){
          fprintf(ferr,"x=%g, y=%g: %g %g %g gh=%g %d\n", x, y, uf.x[-1], uf.x[], uf.x[1], uf.x[ghost], ghost );
          assert(uf.x[ghost] == x);
  }
  fprintf(ferr,"right\n" );

  foreach_boundary(bottom){
          fprintf(ferr,"x=%g, y=%g: %g %g %g\n", x, y, uf.y[-1], uf.y[], uf.y[1] );
          assert(uf.y[] == y);
  }
  fprintf(ferr,"bottom\n" );

  foreach_boundary(top){
          fprintf(ferr,"x=%g, y=%g: %g %g %g\n", x, y, uf.y[-1], uf.y[], uf.y[1] );
          assert(uf.y[ghost] == y);
  }

  foreach_boundary(top){
//      fprintf(ferr, "event2: x=%g y=%g  %g %g %g %g \n", x, y, s[-1], s[], s[1], s[ghost]);
      assert (s[ghost] == y+1);
  }
  return 1;
}

// The setup code:
int main() {
  run();
}
// Hooks for the solver events..
event event1 (i++) {
  uf.n[left] = 22;
  uf.n[right] = 44;
  uf.n[bottom] = 2;
  uf.n[top] = 4;
  s[top] = y;
}

event event2 (i++) {
  uf.n[left] = x;
  uf.n[right] = x;
  uf.n[bottom] = y;
  uf.n[top] = y;
  s[top] = y+1;
}





//#include "run.h"
//
//// The "Solver" code
//face vector uf[];
//uf.n[bottom] = 77;
//uf.n[top] = 99;
//
////void projection(face vector uf){
//////    foreach_face(){
//////        uf.y[]=7;
//////    }
////    boundary({uf});
////}
//
///*event event0 (i++) {
//  boundary({uf});
////  projection(uf);
//  foreach_boundary(bottom){
////    fprintf(ferr, "bottom: %d\n", ghost);
//    assert (uf.y[] == 77);
//  }
//  foreach_boundary(top){
////    fprintf(ferr, "top: %d\n", ghost);
//    fprintf(ferr, "++top: %g %g %g 99\n", uf.y[-1], uf.y[], uf.y[1]);
//    assert (uf.y[1] == 99);
//  }
//}*/
//
//event event1 (i++) {
//boundary({uf});
//// projection(uf);
//foreach_face(){
//    fprintf(ferr,"x=%g, y=%g: %g %g %g\n", x, y, uf.y[-1], uf.y[], uf.y[1] );
//}
//foreach_boundary(bottom)
//assert (uf.y[] == 2);
//fprintf(ferr,"bottom\n" );
//foreach_boundary(top){
//        fprintf(ferr,"top: %g %g %g 4 ", uf.y[-1], uf.y[], uf.y[1] );
//        assert(uf.y[1] == 4);
//}
//
//}
//
//event event2 (i++) {
//boundary({uf});
////projection(uf);
//foreach_boundary(bottom)
//assert (uf.y[] == 1);
//foreach_boundary(top)
//assert (uf.y[1] == 5);
//// return 1;
//}
//
//// The setup code:
//int main() {
//    run();
//}
//// Hooks for the solver events..
//event event1 (i++) {
//uf.n[bottom] = 2;
//uf.n[top] = 4;
//}
//
//event event2 (i++) {
//uf.n[bottom] = 1;
//uf.n[top] = 5;
//}
//
//event event3 (i++) {
////  projection(uf);
////  foreach_boundary(bottom)
////    assert (uf.y[] == 1);
////  foreach_boundary(top)
////    assert (uf.y[1] == 5);
//return 1;
//
//}

