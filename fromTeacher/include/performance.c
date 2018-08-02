#ifndef _PERFORMANCE_C 
#define _PERFORMANCE_C 1 

#include<complex.h>
#include"../../include/complex.c"

extern int N,T,Tp;

int Termination(complex **X_es,complex **X_es_pre)
{
  int n,t;
  double eps;

  eps = 1.0e-2;
  for(n=0;n<N;n++) 
  for(t=0;t<T;t++) {
      if(cabs(X_es[n][t]-X_es_pre[n][t])>eps) return 0;
  }

  return 1;
}


double MSE(complex **H,complex **H_es)
{
  int m,n;
  double mse;

  mse = 0.0;
  for(m=0;m<M;m++) {
    for(n=0;n<N;n++) {
      mse += CAbs2(H[m][n] - H_es[m][n]);
    }
  }
  mse /= (double)N;

  return mse;
}

double BER(int t,complex **X,complex **X_es)
{
  int n;
  double ber;

  ber = 0.0;
  for(n=0;n<N;n++) {
    if(creal(X[n][t])*creal(X_es[n][t])<0.0) ber++;
    if(cimag(X[n][t])*cimag(X_es[n][t])<0.0) ber++;
  }
  ber /= 2.0*N*(T-Tp);

  return ber;
}

#endif
