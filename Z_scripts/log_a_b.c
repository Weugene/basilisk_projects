#include <stdio.h>
#include <math.h>
int main(void)
{
	for (int N = 32; N < 512; N *= 2){
        	int lev = log(N)/log(2);
		printf("log2 N is %d\n", lev);
	}
	return 0;
}


