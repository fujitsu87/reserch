#ifndef _MODULATION_C 
#define _MODULATION_C 1 

#include<complex.h>

#include"../../include/random_number.c"

#define BPSK(x) (1.0 - 2.0*(x))

void QPSKModulation(int k,int t,complex **X)
{
  int i,j;
  double a;

  a = 1.0/sqrt(2.0);
  for(i=0;i<k;i++) {
    for(j=0;j<t;j++) {
      X[i][j] = a*BPSK(UniformBit()) + I*a*BPSK(UniformBit());
      //printf("%g ",cimag(X[i][j]));
    }
    //printf("\n");
  }
}

#endif
