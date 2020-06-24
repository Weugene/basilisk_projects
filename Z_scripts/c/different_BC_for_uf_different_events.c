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
    //face vectors
    foreach_boundary(left){
        assert(uf.x[] == 22);
        fprintf(ferr,"ufx: %g %g %g ufy: %g %g %g \n",uf.x[-1], uf.x[], uf.x[1], uf.y[-1], uf.y[], uf.y[1] );
        assert(uf.y[ghost] == 222);
        assert(uf.y[-1,0] == 222);
    }

    foreach_boundary(right){
        assert(uf.x[ghost] == 44); // they are the same
        assert(uf.x[1,0] == 44);
        assert(uf.y[ghost] == 444);
        assert(uf.y[1,0] == 444);
    }

    foreach_boundary(bottom){
        assert(uf.y[] == 2);
        assert(uf.x[ghost] == 22);
    }

    foreach_boundary(top){
        assert(uf.y[ghost] == 4); // they are the same
        assert(uf.y[0,1] == 4);
        assert(uf.x[0,1] == 44);
    }
    //scalars
    foreach_boundary(bottom){
        assert(s[ghost] == 3); // they are the same
        assert(s[0,-1] == 3);
    }

    foreach_boundary(top){
        assert(s[ghost] == 4); // they are the same
        assert(s[0,1] == 4);
    }

}

event event2 (i++) {
    boundary({s,uf});
    //face vectors
    foreach_boundary(left){
        assert(uf.x[] == x);
    }

    foreach_boundary(right){
        assert(uf.x[ghost] == x);
    }

    foreach_boundary(bottom){
        assert(uf.y[] == y);
    }

    foreach_boundary(top){
        assert(uf.y[ghost] == y);
    }
//scalars
    foreach_boundary(bottom){
        assert(s[ghost] == 1);
        assert(s[0,-1] == 1);
    }

    foreach_boundary(top){
        assert(s[ghost] == 2);
        assert(s[0,1] == 2);
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
    uf.t[left] = 222;
    uf.n[right] = 44;
    uf.t[right] = 444;
    uf.n[bottom] = 2;
    uf.t[bottom] = 22;
    uf.n[top] = 4;
    uf.t[top] = 44;
    s[bottom] = 3;
    s[top] = 4;
}

event event2 (i++) {
    uf.n[left] = x;
    uf.n[right] = x;
    uf.n[bottom] = y;
    uf.n[top] = y;
    s[bottom] = 1;
    s[top] = 2;
}


//#include "run.h"
//
//// The "Solver" code
//face vector uf[];
//scalar s[];
//event event1 (i++) {
//  boundary({uf,s});
////face vectors
//  foreach_boundary(left){
//      assert(uf.x[] == 22);
////      fprintf(ferr,"ufx: %g %g %g ufy: %g %g %g", uf.x[-1], uf.x[], uf.x[1], uf.y[-1], uf.y[], uf.y[1]);
////      assert(uf.y[-1] == 222);
//  }
//
//  foreach_boundary(right){
//      assert(uf.x[ghost] == 44); // they are the same
//      assert(uf.x[1,0] == 44);
//  }
//
//  foreach_boundary(bottom){
//      assert(uf.y[] == 2);
//  }
//
//  foreach_boundary(top){
//      assert(uf.y[ghost] == 4); // they are the same
//      assert(uf.y[0,1] == 4);
//  }
////scalars
//  foreach_boundary(bottom){
//      assert(s[ghost] == y+3); // they are the same
//      assert(s[0,-1] == y+3);
//  }
//
//  foreach_boundary(top){
//      assert(s[ghost] == y); // they are the same
//      assert(s[0,1] == y+3);
//  }
//
//}
//
//event event2 (i++) {
//  boundary({s,uf});
////face vectors
//  foreach_boundary(left){
//      assert(uf.x[] == x);
//  }
//
//  foreach_boundary(right){
//      assert(uf.x[ghost] == x);
//  }
//
//  foreach_boundary(bottom){
//      assert(uf.y[] == y);
//  }
//
//  foreach_boundary(top){
//      assert(uf.y[ghost] == y);
//  }
////scalars
//  foreach_boundary(bottom){
//      assert(s[ghost] == 1);
//      assert(s[0,-1] == 1);
//  }
//
//  foreach_boundary(top){
//      assert(s[ghost] == 2);
//      assert(s[0,1] == 2);
//  }
//  return 1;
//}
//
//// The setup code:
//int main() {
//  run();
//}
//// Hooks for the solver events..
//event event1 (i++) {
//  uf.n[left] = 22;
////  uf.t[left] = 222;
//  uf.n[right] = 44;
//  uf.n[bottom] = 2;
//  uf.n[top] = 4;
//  s[bottom] = 3;
//  s[top] = 4;
//}
//
//event event2 (i++) {
//  uf.n[left] = x;
//  uf.n[right] = x;
//  uf.n[bottom] = y;
//  uf.n[top] = y;
//  s[bottom] = 1;
//  s[top] = 2;
//}




//#include "run.h"
//
//// The "Solver" code
//face vector uf[];
//scalar s[];
//event event1 (i++) {
//boundary({uf,s});
//
//foreach_boundary(left){
////          fprintf(ferr,"x=%g, y=%g: %g %g %g gh=%g %d\n", x, y, uf.x[-1], uf.x[], uf.x[1], uf.x[ghost], ghost );
//        assert(uf.x[] == 22);
//}
//fprintf(ferr,"left\n" );
//
//foreach_boundary(right){
////          fprintf(ferr,"x=%g, y=%g: %g %g %g gh=%g %d\n", x, y, uf.x[-1], uf.x[], uf.x[1], uf.x[ghost], ghost );
//        assert(uf.x[1] == 44);
//}
//fprintf(ferr,"right\n" );
//
//foreach_boundary(bottom){
////        fprintf(ferr,"x=%g, y=%g: %g %g %g\n", x, y, uf.y[-1], uf.y[], uf.y[1] );
//        assert(uf.y[0,0] == 2);
//}
//fprintf(ferr,"bottom\n" );
//
//foreach_boundary(top){
////      fprintf(ferr,"x=%g, y=%g: %g %g %g\n", x, y, uf.y[0,-1], uf.y[], uf.y[0,1] );
//        assert(uf.y[0,1] == 4);
//}
//fprintf(ferr,"top\n" );
//
//foreach_boundary(top){
////    fprintf(ferr, "event1: x=%g y=%g | %g %g %g %g \n", x, y, s[-1], s[], s[1], s[ghost]);
//        fprintf(ferr, "event1: ig=%d gh=%d %d %g %g\n", ig, ghost, s[ghost], s[0]);
//        assert(s[0,1] == y);
//}
//
//}
//
//event event2 (i++) {
//boundary({s,uf});
//foreach_boundary(left){
//        fprintf(ferr,"x=%g, y=%g: %g %g %g gh=%g %d\n", x, y, uf.x[-1], uf.x[], uf.x[1], uf.x[ghost], ghost );
//        assert(uf.x[] == x);
//}
//fprintf(ferr,"left\n" );
//
//foreach_boundary(right){
//        fprintf(ferr,"x=%g, y=%g: %g %g %g gh=%g %d\n", x, y, uf.x[-1], uf.x[], uf.x[1], uf.x[ghost], ghost );
//        assert(uf.x[1] == x);
//}
//fprintf(ferr,"right\n" );
//
//foreach_boundary(bottom){
//        fprintf(ferr,"x=%g, y=%g: %g %g %g\n", x, y, uf.y[-1], uf.y[], uf.y[1] );
//        assert(uf.y[0,0] == y);
//}
//fprintf(ferr,"bottom\n" );
//
//foreach_boundary(top){
//        fprintf(ferr,"x=%g, y=%g: %g %g %g\n", x, y, uf.y[-1], uf.y[], uf.y[1] );
//        assert(uf.y[0,1] == y);
//}
//
//foreach_boundary(top){
////      fprintf(ferr, "event2: x=%g y=%g  %g %g %g %g \n", x, y, s[-1], s[], s[1], s[ghost]);
//        assert (s[0,1] == y+1);
//}
//return 1;
//}




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

