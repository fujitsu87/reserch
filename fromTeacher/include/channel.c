#ifndef _CHANNEL_C 
#define _CHANNEL_C 1 

#include<complex.h>
#include<math.h>

#include"../../include/random_number.c"

extern int M,N,T;
extern double sigma;

void ComplexGaussianMatrix(int m,int n,complex **H,double s)
{
  int i,j;

  for(i=0;i<m;i++) {
    for(j=0;j<n;j++) {
      H[i][j] = ComplexGaussianNoise(s);
      //printf("%g ",creal(H[i][j]));
    }
    //printf("\n");
  }
}

void MIMOChannel(complex **Y,complex **H,complex **X)
{
  int m,n,t;

  for(t=0;t<T;t++) {
    for(m=0;m<M;m++) {
      Y[m][t] = 0.0;
      for(n=0;n<N;n++) {
	Y[m][t] += H[m][n]*X[n][t]; 
      }
      Y[m][t] += ComplexGaussianNoise(sigma);
      //printf("%g ",creal(Y[m][t]));
    }
    //printf("\n");
  }
}


#endif