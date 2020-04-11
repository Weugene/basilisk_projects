#include "run.h"
vector b[];
(const) vector bU ;
int main(int argc, char * argv[])
{
    init_grid(16);
    const vector abc[] = {1.,0.,0.};
    bU = abc;

//    foreach_dimension(){
//        foreach() {
//            b.x[] = bU.x[];
//        }
//    }
  
    foreach() {
        foreach_dimension(){
            b.x[] = bU.x[];
        }
    }
    foreach() fprintf(ferr,"b=%g %g bU=%g %g\n", b.x[], b.y[], bU.x[], bU.y[]);

return 0;
}
