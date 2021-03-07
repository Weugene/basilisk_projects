#include <stdio.h>
#include <math.h>
int main(void)
{

for (double i=-1; i<=1; i++)
    for(double j=-1; j<=1;j++){
        printf("atan2 of x=%f y=%f is %f\n", i, j, 180*atan2(j, i)/pi);
}
return 0;
}
