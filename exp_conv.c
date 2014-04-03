#include <math.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>


int main(int argc, char *argv[]){

/* Numerical parameters - default values */

int N = (int)pow(2,22); /* Total number of grid points; so 1D arrays will run [0:Nx-1] */
int Nt = 1; /* Number of time steps */


/* Check inputs */
if (argc>1){
	int base = atoi(argv[1]);
	int power = atoi(argv[2]);
	N = (int)pow(base,power);
}

if (argc>3){
	Nt = atoi(argv[3]);
}


printf("Number of grid points: %d\n", N);
printf("Number of time steps: %d\n",Nt);

/* Domain set up */
float L = 1.0; /* Total length of the domain */
float dx = L/((float)(N-1));
float nu = 1.0;
float alpha = nu/dx;
float ex = exp(-nu);
float eL = exp(-alpha*L);

/* Quadrature weights */
float P = 1.0 - (1.0-ex)/nu;
float Q = -ex+(1.0-ex)/nu;
float R = (1.0-ex)/(nu*nu)-(1.0+ex)/(2.0*nu);





/* Solution array set up */

/* Main solution arrays */
float integrand[N];
float I[N];


float JR[N],JL[N],IR[N],IL[N];


/* Some indices for utility */
int n;
int j;
int i;


for(j=0; j<N; j++){
	integrand[j]=1.0;
	IL[j]=0.0;
	IR[j]=0.0;
	JL[j]=0.0;
	JR[j]=0.0;
	I[j]=0.0;
}



/*************************************************** Time loop **************************************************************************/
clock_t start = clock(), diff;
for(n = 1; n <= Nt; n++){

	for(j=1; j<N-1; j++){
		JL[j]=P*integrand[j]+Q*integrand[j-1]+R*(integrand[j+1]-2.0*integrand[j]+integrand[j-1]);
	}
	
	JL[N-1]=P*integrand[N-1]+Q*integrand[N-2]+R*(integrand[1]-2.0*integrand[N-1]+integrand[N-2]);	
	JL[0]=(float) 0;


	for(j=0; j<N-1; j++){
		IL[j+1]=JL[j+1]+ex*IL[j];
	}

	for(j=1; j<N-1; j++){
		JR[j]=P*integrand[j]+Q*integrand[j+1]+R*(integrand[j+1]-2.0*integrand[j]+integrand[j-1]);
	}
	
	JR[0]=P*integrand[0]+Q*integrand[1]+R*(integrand[1]-2.0*integrand[0]+integrand[N-2]);	
	JR[N-1]=(float) 0;


	for(j=1; j<N; j++){
		IR[N-j-1]=JR[N-j-1]+ex*IR[N-j];
	}
	
	
	for(j=0; j<N; j++){
		I[j]=IL[j]+IL[N-1-j];
	}
}

diff = clock() - start;
int sec = diff/ CLOCKS_PER_SEC;
printf("Time taken: %d seconds\n", sec);

/* Compute test integral value. */

float x;
float err = 0.0;
float err_temp;
for(j=0; j<N; j++){
	x = (float)j*dx;
	err_temp = abs(I[j]-(2.0-exp(-alpha*x)-exp(-alpha*(L-x))));
	
	if (err_temp>err)
	{
		err = err_temp;
	}
	
}


printf("Maximum error: %d \n", err);

return 0;
}
