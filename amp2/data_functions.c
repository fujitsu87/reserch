

double Computation_Vx(int k,int t)
{
  int n;
  double vx;

  vx = 0.0;
  for(n=0;n<N;n++) {
    vx += gsl_complex_abs2(gsl_matrix_complex_get(h_h,n,k))/gsl_matrix_get(zeta_d,n,t);
  }
  vx = 1.0/vx;

  return vx;
}


gsl_complex Computation_X_es(int k,int t,gsl_matrix_complex *X_es_pre)
{
  int n;
  gsl_complex a,b,ans;

  a.dat[0] = 0.0;
  a.dat[1] = 0.0; 
  for(n=0;n<N;n++) {
    // a += conj(H_es[m][n])*Z_d[m][t]/zeta_d[m][t];
    a = gsl_complex_add(
        a,
        gsl_complex_div_real(
            gsl_complex_mul(
                gsl_complex_conjugate(
                    gsl_matrix_complex_get(h_h,n,k)
                ),
                gsl_matrix_complex_get(z_d,n,t)
            ),
            gsl_matrix_get(zeta_d,n,t)
        )
    );
  }
  a = gsl_complex_mul_real(
    a,
    gsl_matrix_get(xi,k,t)
  );

  b.dat[0] = 0.0;
  b.dat[1] = 0.0; 
  for(n=0;n<N;n++) {
    // b += Vh[m][n]*conj(Z_c[m][t])*Z_d[m][t]/(zeta_c[m][t]*zeta_d[m][t]);
    b = gsl_complex_add(
        b,
        gsl_complex_div_real(
            gsl_complex_mul_real(
                gsl_complex_mul(
                    gsl_complex_conjugate(
                        gsl_matrix_complex_get(z_c,n,t)
                    ),
                    gsl_matrix_complex_get(z_d,n,t)
                ),
                gsl_matrix_get(eta,n,k)
            ),
            gsl_matrix_get(zeta_d,n,t)*gsl_matrix_get(zeta_c,n,t)
        )
    );
  }
//   b *= Vx[n][t]*X_es[n][t];
  b = gsl_complex_mul(
      b,
      gsl_complex_mul_real(
          gsl_matrix_complex_get(x_h,k,t),
          gsl_matrix_get(xi,k,t)
      )
  );
    ans = gsl_complex_add(
        gsl_matrix_complex_get(x_h,k,t),
        a
    );
    ans = gsl_complex_sub(
        ans,
        b
  );
  return ans;
}
gsl_complex SoftDecision(gsl_complex x_es,double v)
{
  double a,b;
  gsl_complex ans;
  a = 1.0/sqrt(2.0); b = 2.0*a/v;

    //a*tanh(b*creal(x_es)) + I*a*tanh(b*cimag(x_es));
  ans.dat[0] = a*tanh(b*x_es.dat[0]);
  ans.dat[1] = a*tanh(b*x_es.dat[1]);
 
  return ans;
}
//------------------データ推定器----------------------------------
void data_estimation(FILE *fp_x,int t,double *ber)
{
    int n,k,i;
    double a;
    // complex **X_es_pre,**H_es_pre,*x_es_pre;

    gsl_matrix_complex *X_es_pre,*H_es_pre;
    gsl_vector_complex *x_es_pre;

    X_es_pre = gsl_matrix_complex_calloc(K,T);
    x_es_pre = gsl_vector_complex_calloc(K);
    H_es_pre = gsl_matrix_complex_calloc(N,K);

    for(k=0;k<K;k++) {
        gsl_matrix_complex_set(
            X_es_pre,k,t,gsl_matrix_complex_get(x_h,k,t)
        );
        gsl_vector_complex_set(
            x_es_pre,k,gsl_matrix_complex_get(x_h,k,t)
        ); 
    }
    gsl_matrix_complex_memcpy(H_es_pre,h_h);
    for(i=0;i<X_LOOP;i++) {
        for(n=0;n<N;n++) {
            gsl_complex z = Computation_z(n,t,X_es_pre,H_es_pre);
            gsl_matrix_complex_set(z_d,n,t,z);
            double zeta = Computation_zeta(n,t); 
            gsl_matrix_set(zeta_d,n,t,zeta);
        }
        for(k=0;k<K;k++) {
            if(gsl_matrix_get(pilot,k,t) == 0){
                gsl_matrix_set(
                    xi,k,t,
                    Computation_Vx(k,t)
                );
                gsl_matrix_complex_set(
                    x_h,k,t,
                    Computation_X_es(k,t,X_es_pre)
                );
            }
        }
        for(k=0;k<K;k++) {
            if(gsl_matrix_get(pilot,k,t) == 0){
                gsl_matrix_complex_set(
                    x_h,k,t,
                    SoftDecision(gsl_matrix_complex_get(x_h,k,t),gsl_matrix_get(xi,k,t)
                    )
                );
                gsl_matrix_set(
                    xi,k,t,1.0 - gsl_complex_abs2(gsl_matrix_complex_get(x_h,k,t))
                );
            }
        }
        a = BER(t);
        ber[i] += a;
        for(k=0;k<K;k++){
            if(gsl_matrix_get(pilot,k,t) == 0){
                gsl_vector_complex_set(
                    x_es_pre,
                    k,
                    gsl_matrix_complex_get(x_h,k,t)
                );
            }
        }  
    }
    X_es_pre = GSLMatrixFree(X_es_pre);
    H_es_pre = GSLMatrixFree(H_es_pre);
    x_es_pre = GSLVectorFree(x_es_pre);
    // x_b = GSLMatrixFree(x_b);
    // PrintMatrix(stdout,K/2,T,x);
    // PrintMatrix(stdout,K/2,T,x_h);
}