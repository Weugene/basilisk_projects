#include "run.h"
int main()
{
	size (1.0);
	origin (-0.5*L0, -0.5*L0);
	periodic (right);
	periodic (top);
	N = 1 << 4;
	run();
}
event printtt(i++; i<1){
fprintf(ferr, "X0=%g Y0=%g L0=%g", X0,Y0,L0);
foreach() fprintf(ferr, "cell: x=%.15g y=%.15g Delta=%.15g \n", x, y, Delta);
foreach_vertex() fprintf(ferr, "vert: x=%.15g y=%.15g Delta=%.15g \n ", x, y, Delta);


}
