#ifndef _AMP_OPERATION_C 
#define _AMP_OPERATION_C 1 

#include<complex.h>
#include<gsl/gsl_cdf.h>

#include"../../include/vector.c"
#include"../../include/complex.c"

extern int M,N,T;
extern double sigma2;

double Computation_zeta(int m,int t,complex **H_es,double **Vh,
			complex **X_es,double **Vx)
{
  int n;
  double zeta;

  zeta = sigma2;
  for(n=0;n<N;n++) {
    zeta += Vh[m][n]*(Vx[n][t] + CAbs2(X_es[n][t]))
      + CAbs2(H_es[m][n])*Vx[n][t];
  }

  return zeta;
}

complex Computation_z(int m,int t,complex **Y,complex **H_es,double **Vh,
		      complex **X_es,double **Vx,complex **Z_c,double **zeta_c,
		      complex **Z_d,double **zeta_d,complex **X_es_pre,
		      complex **H_es_pre)
{
  int n;
  complex a,b,c;

  a = 0.0;
  for(n=0;n<N;n++) {
    a += H_es[m][n]*X_es[n][t];
  }

  b = 0.0;
  for(n=0;n<N;n++) {
    b += Vh[m][n]*conj(X_es_pre[n][t])*X_es[n][t];
  }

  c = 0.0;
  for(n=0;n<N;n++) {
    c += conj(H_es_pre[m][n])*H_es[m][n]*Vx[n][t];
  }

  //printf("  %g %g\n",creal(Z_c[m][t])/zeta_c[m][t],zeta_c[m][t]);

  return Y[m][t] - a + b*Z_c[m][t]/zeta_c[m][t] 
    + c*Z_d[m][t]/zeta_d[m][t]; 
}




#endif
