#include "stdio.h"
#include "math.h"

//David Schmidt

int main (void)
{
    float x[5] = {0, 0, 0.6, 1.2, 1.8};  // x[1..4] used
    float f[5] = {0, 0, 1.0, 2.0, 5.0};  // f[1..4] used
    float u[5];                          // u[1..4] used
    float c[5];                          // c[1..4] used

    int i,m;     //loop indices
    int N = 4;   // number of points


    for (i=1; i< N+1; i++)
    {
	u[i]=f[i];
    }

    c[1] = u[1];
    printf(" Coefficient %d is %f \n",1,c[1]);

    for (m = 2; m<N+1; m++) // run over columns
    {
	for (i = 1; i<=N-m+2; i++) // run over rows
	{
	    u[i] = (u[i+1]-u[i])/(x[i+m-1]-x[i]);
	}

	c[m] = u[1];
	printf(" Coefficient %d is %f \n",m,c[m]);
    }

    return(0);
}
