#ifndef _DATA_ESTIMATION_C 
#define _DATA_ESTIMATION_C 1 

#include<math.h>
#include<complex.h>

#include"../../include/complex.c"
#include"./AMP_operation.c"

extern int M,N,T;

int Termination_d(int t,complex **X_es,complex *x_es_pre)
{
  int n;
  double eps;

  eps = 1.0e-2;
  for(n=0;n<N;n++) {
      if(cabs(X_es[n][t]-x_es_pre[n])/cabs(x_es_pre[n])>eps) return 0;
  }

  return 1;
}

double Computation_Vx(int n,int t,complex **H_es,double **zeta_d)
{
  int m;
  double vx;

  vx = 0.0;
  for(m=0;m<M;m++) {
    vx += CAbs2(H_es[m][n])/zeta_d[m][t];
    //printf("  %g %g\n",M*CAbs2(H_es[m][n]),zeta_d[m][t]);
  }
  vx = 1.0/vx;

  return vx;
}


complex Computation_X_es(int n,int t,complex **H_es,double **Vh,
			 complex **X_es,double **Vx,complex **Z_d,
			 double **zeta_d,complex **Z_c,double **zeta_c,
			 complex **X_es_pre)
{
  int m;
  complex a,b;

  a = 0.0; 
  for(m=0;m<M;m++) {
    a += conj(H_es[m][n])*Z_d[m][t]/zeta_d[m][t];
  }
  a *= Vx[n][t];

  b = 0.0; 
  for(m=0;m<M;m++) {
    b += Vh[m][n]*conj(Z_c[m][t])*Z_d[m][t]/(zeta_c[m][t]*zeta_d[m][t]);
  }
  b *= Vx[n][t]*X_es[n][t];

  return X_es[n][t] + a - b;
}

complex SoftDecision(complex x_es,double v)
{
  double a,b;

  a = 1.0/sqrt(2.0); b = 2.0*a/v;

  return a*tanh(b*creal(x_es)) + I*a*tanh(b*cimag(x_es));
}

void Data_Estimation(int t,int iter,complex **Y,complex **X_es,double **Vx,
		     complex **Z_d,double **zeta_d,complex **H_es,double **Vh,
		     complex **Z_c,double **zeta_c,complex **X,double *ber)
{
  int m,n,i;
  double a;
  complex **X_es_pre,**H_es_pre,*x_es_pre;

  X_es_pre = CPMatrixAlloc(0,N-1,0,T-1);
  x_es_pre = CPVectorAlloc(0,N-1);
  H_es_pre = CPMatrixAlloc(0,M-1,0,N-1);

  for(n=0;n<N;n++) {
    X_es_pre[n][t] = x_es_pre[n] = X_es[n][t];
  }

  for(m=0;m<M;m++) {
    for(n=0;n<N;n++) {
      H_es_pre[m][n] = H_es[m][n];
    }
  }

  for(i=0;i<iter;i++) {
    for(m=0;m<M;m++) {
      Z_d[m][t] = Computation_z(m,t,Y,H_es,Vh,X_es,Vx,Z_c,zeta_c,Z_d,zeta_d,
				X_es_pre,H_es_pre);
      //printf("%g %g\n",creal(Z_d[m][t]),cimag(Z_d[m][t]));
      zeta_d[m][t] = Computation_zeta(m,t,H_es,Vh,X_es,Vx); 
    }

    for(n=0;n<N;n++) {
      Vx[n][t] = Computation_Vx(n,t,H_es,zeta_d);
      X_es[n][t] = Computation_X_es(n,t,H_es,Vh,X_es,Vx,
				    Z_d,zeta_d,Z_c,zeta_c,X_es_pre);
      //printf("%g %g %g\n",creal(X_es[n][t]),cimag(X_es[n][t]),Vx[n][t]);
    }

    for(n=0;n<N;n++) {
      X_es[n][t] = SoftDecision(X_es[n][t],Vx[n][t]);
      Vx[n][t] = 1.0 - CAbs2(X_es[n][t]);
    }

    a = BER(t,X,X_es);
    if(Termination_d(t,X_es,x_es_pre)!=0) {
      for(;i<iter;i++) {
        ber[i] = +a;
      }
      break;
    }
    else {
      ber[i] = +a;
      for(n=0;n<N;n++) x_es_pre[n] = X_es[n][t]; 
    }

  }

  X_es_pre = CPMatrixFree(X_es_pre,0,N-1,0,T-1);
  x_es_pre = CPVectorFree(x_es_pre,0,N-1);
  H_es_pre = CPMatrixFree(H_es_pre,0,M-1,0,N-1);
}


#endif
