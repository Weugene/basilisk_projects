#include "run.h"
vector b[];
const vector zerocf[] = {1.,0.,0.};
int main(int argc, char * argv[])
{
    init_grid(16);
    (const) vector bU = zerocf; 
    foreach_dimension(){
        foreach() {
            b.x[] = bU.x[];
        }
    }
    foreach() fprintf(ferr,"b=%g %g bU=%g %g\n", b.x[], b.y[], bU.x[], bU.y[]);

return 0;
}
