#include "run.h"

// The "Solver" code
scalar s[];

event event1 (i++) {
//foreach_level_or_leaf(2) {
//        s[] = x;
//        fprintf(ferr, "%g\n", s[]);
////        return 1;
//    }

    double pcor;

    int k=0;
    foreach(){
        if(k==0){
            pcor = s[];
//            k++;
            break;
        }
    }
    fprintf(ferr, "event1 pcor=%g", pcor);
//    fprintf(ferr, "%g\n", ((double *) ((((Tree *)grid)->L[point.level]->m[point.i+0][point.j+0]) + sizeof(Cell)))[(s.i)]);
    return 1;
}


// The setup code:
int main() {
    run();
}

