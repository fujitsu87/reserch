#ifndef _CHANNEL_ESTIMATION_C 
#define _CHANNEL_ESTIMATION_C 1 

#include<complex.h>

#include"../../include/complex.c"
#include"../../include/vector.c"

#include"./AMP_operation.c"
#include"./performance.c"

extern int M,N,T;
extern double sigma2;

int Termination_c(complex **H_es,complex **H_es_pre)
{
  int m,n;
  double eps;

  eps = 1.0e-2;
  for(m=0;m<M;m++) {
    for(n=0;n<N;n++) {
      if(cabs(H_es[m][n]-H_es_pre[m][n])/cabs(H_es_pre[m][n])>eps) return 0;
    }
  }

  return 1;
}

double Computation_Vh(int m,int n,complex **X_es,double **zeta_c)
{
  int t;
  double vh;

  vh = 0.0;
  for(t=0;t<T;t++) {
    vh += CAbs2(X_es[n][t])/zeta_c[m][t];
  }
  vh /= (double)M;
  vh = 1.0/(M*(1.0 + vh));

  return vh;
}

complex Computation_H_es(int m,int n,complex **X_es,double **Vx,
			 complex **H_es,double **Vh,complex **Z_c,
			 double **zeta_c,complex **Z_d,double **zeta_d,
			 complex **H_es_pre)
{
  int t;
  complex a,b;

  a = 0.0; 
  for(t=0;t<T;t++) {
    a += conj(X_es[n][t])*Z_c[m][t]/zeta_c[m][t];
    //printf("  %g\n",creal(Z_c[m][t]/zeta_c[m][t]));
  } 
  a *= Vh[m][n];

  b = 0.0; 
  for(t=0;t<T;t++) {
    b += Vx[n][t]*conj(Z_d[m][t])*Z_c[m][t]/(zeta_d[m][t]*zeta_c[m][t]);
  }
  b *= Vh[m][n]*H_es_pre[m][n];

  return (1.0-M*Vh[m][n])*H_es[m][n] + a - b;
}


void Channel_Estimation(int iter,complex **Y,complex **H_es,double **Vh,
			complex **Z_c,double **zeta_c,complex **X_es,
			double **Vx,complex **Z_d,double **zeta_d,
			complex **H,double *mse)
{
  int m,n,t,i;
  double a;
  complex **X_es_pre,**H_es0,**H_es_pre;

  X_es_pre = CPMatrixAlloc(0,N-1,0,T-1);
  H_es0 = CPMatrixAlloc(0,M-1,0,N-1);
  H_es_pre = CPMatrixAlloc(0,M-1,0,N-1);

  for(n=0;n<N;n++) {
    for(t=0;t<T;t++) X_es_pre[n][t] = X_es[n][t];
  }

  for(m=0;m<M;m++) {
    for(n=0;n<N;n++) {
      H_es0[m][n] = H_es_pre[m][n] = H_es[m][n];
    }
  }

  for(i=0;i<iter;i++) {
    for(m=0;m<M;m++) {
      for(t=0;t<T;t++) {
	Z_c[m][t] = Computation_z(m,t,Y,H_es,Vh,X_es,Vx,Z_c,zeta_c,Z_d,zeta_d,X_es_pre,H_es0);
	zeta_c[m][t] = Computation_zeta(m,t,H_es,Vh,X_es,Vx);
      }
      //printf("%d %g %g\n",m,creal(Z_c[m][t]),zeta_c[m][t]-sigma2);
    }
    //printf("\n");


    for(m=0;m<M;m++) {
      for(n=0;n<N;n++) {
	Vh[m][n] = Computation_Vh(m,n,X_es,zeta_c);
	H_es[m][n] = Computation_H_es(m,n,X_es,Vx,H_es,Vh,
				    Z_c,zeta_c,Z_d,zeta_d,H_es0);
	//printf("%g %g %g %g %g\n",sqrt(M)*creal(H[m][n]),sqrt(M)*creal(H_es[m][n]),sqrt(M)*cimag(H[m][n]),sqrt(M)*cimag(H_es[m][n]),M*Vh[m][n]);
      }
    }
    //printf("\n");
    
    a = MSE(H,H_es);
    if(Termination_c(H_es,H_es_pre)!=0) {
      //printf("%d\n",i);
      for(;i<iter;i++) mse[i] += a;
      break;
    }
    else {
      mse[i] += a;
      for(m=0;m<M;m++) {
	for(n=0;n<N;n++) H_es_pre[m][n] = H_es[m][n];
      } 
    }
  }

  X_es_pre = CPMatrixFree(X_es_pre,0,N-1,0,T-1);
  H_es0 = CPMatrixFree(H_es0,0,M-1,0,N-1);
  H_es_pre = CPMatrixFree(H_es_pre,0,M-1,0,N-1);
}

#endif
