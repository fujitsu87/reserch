//------------------culc_eta------------------------------------
void culc_eta()
{
    int n,k,t;
    double tmp1,tmp2;

    for (n = 0; n < N; ++n)
    {
        for (k = 0; k < K; ++k)
        {
            tmp2 = 0;
            for (t = 0; t < T; ++t)
            {
                tmp2 += gsl_matrix_get(x_h_abs2,k,t)/(N0+ gsl_matrix_get(zeta,n,t));
            }
            tmp1 = 1 + tmp2/(double)N;

            gsl_matrix_set(eta,n,k,1/tmp1);
        }
    }
    // PrintRealMatrix(stdout,N,K,eta);
}
//------------------culc_h_h-----------------------------------
void culc_h_h()
{

    int n,t,k;
    gsl_complex tmp1;
    gsl_complex tmp2;
    gsl_complex tmp3;
    gsl_complex sec1;
    gsl_complex sec2;
    gsl_complex sec3;
    GSL_SET_COMPLEX(&sec1,0,0);
    GSL_SET_COMPLEX(&sec2,0,0);
    GSL_SET_COMPLEX(&sec3,0,0);

    gsl_matrix_complex *x_h_conju;
    x_h_conju = gsl_matrix_complex_calloc(T,K);
    ConjugateTranspose(x_h,x_h_conju);

    gsl_matrix_complex *h_h_conju;
    h_h_conju = gsl_matrix_complex_calloc(K,N);
    ConjugateTranspose(h_h,h_h_conju);

    double n_sqrt = sqrt((double)N);
    double tmp4,tmp5;
    for (n = 0; n < N; ++n)
    {
        for (k = 0; k < K; ++k)
        {
            a_k = sqrt(user_p[k]/(2*rho_k));
            GSL_SET_COMPLEX(&sec1,0,0);
            GSL_SET_COMPLEX(&sec2,0,0);
            GSL_SET_COMPLEX(&sec3,0,0);
            tmp4 = 0;
            for (t = 0; t < T; ++t)
            {
                GSL_SET_COMPLEX(&tmp1,0,0);
                
                //z * x_h_conju　第1項の計算
                tmp1 = gsl_complex_mul(
                            gsl_matrix_complex_get(x_h_conju,t,k)
                            ,gsl_matrix_complex_get(z,n,t)
                        );
                sec1.dat[0] += tmp1.dat[0];
                sec1.dat[1] += tmp1.dat[1];
            
                //第3項の計算
                GSL_SET_COMPLEX(&tmp2,0,0);
                
                //Aktの計算
                /*
                tmp2 = A(
                            gsl_complex_mul(
                                gsl_matrix_complex_get(h_h_conju,k,n),
                                gsl_matrix_complex_get(z,n,t)
                            )
                            ,k
                            ,t
                        );
                tmp2 = gsl_complex_conjugate(tmp2);
                tmp2 = gsl_complex_mul_real(
                            tmp2
                            ,gsl_matrix_get(xi_b,k,t)
                        );
                tmp2 = gsl_complex_mul(
                            gsl_matrix_complex_get(z,n,t)
                            ,tmp2
                        );
                sec3.dat[0] += tmp2.dat[0];
                sec3.dat[1] += tmp2.dat[1];
                */

                // xi*|z|^2
                tmp4 = tmp4 + gsl_matrix_get(xi,k,t) 
                        *gsl_complex_abs2(gsl_matrix_complex_get(z,n,t));
            }
            sec1.dat[0] = sec1.dat[0]/n_sqrt;
            sec1.dat[1] = sec1.dat[1]/n_sqrt;
            
            sec1 = gsl_complex_mul_real(sec1,gsl_matrix_get(eta,n,k));

            sec2 = gsl_complex_mul_real(
                    gsl_matrix_complex_get(h_h,n,k)
                    ,(1-gsl_matrix_get(eta,n,k))
                );

            // sum(xi*|z|^2)*h_h
            sec3 = gsl_complex_mul_real(
                            gsl_matrix_complex_get(h_h,n,k)
                            ,tmp4
                        );
            sec3.dat[0] = sec3.dat[0]/(double)N;
            sec3.dat[1] = sec3.dat[1]/(double)N;
            sec3 = gsl_complex_mul_real(sec3,(-1)*gsl_matrix_get(eta,n,k));

            GSL_SET_COMPLEX(&tmp3,0,0);
            tmp3 = gsl_complex_add(sec1,sec2);
            tmp3 = gsl_complex_add(tmp3,sec3);

            tmp3.dat[0] = (1.0 - a_h)*GSL_REAL(gsl_matrix_complex_get(h_h,n,k))
                        +a_h*tmp3.dat[0];
            tmp3.dat[1] = (1.0 - a_h)*GSL_IMAG(gsl_matrix_complex_get(h_h,n,k))
                        +a_h*tmp3.dat[1];

            gsl_matrix_complex_set(h_h,n,k,tmp3);
        }
    }
    x_h_conju = GSLMatrixFree(x_h_conju);
    h_h_conju = GSLMatrixFree(h_h_conju);
}

//------------------通信路推定器----------------------------------
void channel_estimation(FILE *fp_x)
{
    culc_zata();
    culc_I();
    culc_z();

    culc_eta();
    culc_h_h();

    // gsl_matrix_set_zero(eta);
    // gsl_matrix_complex_memcpy(h_h,h);
 
    // PrintMatrix(stdout,N,K,h_h);
}
