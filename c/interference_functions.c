void abs2_matrix(gsl_matrix_complex *x, gsl_matrix *y,int s1,int s2);

//-----------zetaの計算-----------------------------------

void culc_zata()
{
    int n,t,k;

    abs2_matrix(x_h,x_h_abs2,K,T);
    abs2_matrix(h_h,h_h_abs2,N,K);

    double sec1,sec2,sec3;

    // tmp ... sum( eta*xi + |h_h|^2*xi eta*|x_h|^2)
    double tmp = 0;
    for (n = 0; n < N; ++n)
    {
        for (t = 0; t < T; ++t)
        {
            tmp = 0;
            for (k = 0; k < K; ++k)
            {
                sec1 = gsl_matrix_get(eta,n,k) * gsl_matrix_get(xi,k,t);
                sec2 = gsl_matrix_get(eta,n,k) * gsl_matrix_get(x_h_abs2,k,t);
                sec3 = gsl_matrix_get(h_h_abs2,n,k) * gsl_matrix_get(xi,k,t);
                tmp += sec1 + sec2 + sec3;
            }
            tmp = tmp/(double)N;
            if(tmp < 0)printf("%f\n",tmp );
            gsl_matrix_set(zeta,n,t,tmp);
        }
    }

}

//----------------I_bの計算--------------------------------
void culc_I()
{

    abs2_matrix(x_h,x_h_abs2,K,T);
    abs2_matrix(h_h,h_h_abs2,N,K);

    int n,t,k;
    // tmp ... 1/n_sqrt + i0
    gsl_complex tmp;

    //tmp2 ... sum ((h_h^2*xi*z) + (eta*x_h^2*z))
    gsl_complex tmp2;

    //tmp3 ... h_h^2*xi*z tmp4 ... eta*x_h^2*z
    gsl_complex tmp3;
    gsl_complex tmp4;
    
    abs2_matrix(x_h,x_h_abs2,K,T);
    abs2_matrix(h_h,h_h_abs2,N,K);

    gsl_matrix_complex *h_h_conju;
    h_h_conju = gsl_matrix_complex_calloc(K,N);
    ConjugateTranspose(h_h,h_h_conju);

    // 1/sqrt(N) h_h * x_h
    AB(h_h,x_h,I_b);
    double n_sqrt = sqrt((double)N);
    GSL_SET_COMPLEX(&tmp,1/n_sqrt,0);
    gsl_matrix_complex_scale(I_b,tmp);

    //xi_b * A()
    for (n = 0; n < N; ++n)
    {
        for (t = 0; t < T; ++t)
        {
            GSL_SET_COMPLEX(&tmp2,0,0);
            for (k = 0; k < K; ++k)
            {
                a_k = sqrt(user_p[k]/(2*rho_k));
                // GSL_SET_COMPLEX(&tmp3,0,0);
               
               /* 旧ver
                //xiを取り出す
                GSL_SET_COMPLEX(&tmp3,gsl_matrix_get(xi_b,k,t),0);
                //Aktの計算
                tmp4 = A(
                            gsl_complex_mul(
                                gsl_matrix_complex_get(h_h_conju,k,n),
                                gsl_matrix_complex_get(z,n,t)
                            )
                            ,k
                            ,t
                        );

                //xi_b * A
                tmp3 = gsl_complex_mul(tmp3,tmp4);
                //h_hとかけ合わせる
                tmp3 = gsl_complex_mul(gsl_matrix_complex_get(h_h,n,k),tmp3);
            */
                //h_h^2 * xi * znt
                tmp3 = gsl_complex_mul_real(
                        gsl_matrix_complex_get(z,n,t)
                        ,gsl_matrix_get(xi,k,t)*gsl_matrix_get(h_h_abs2,n,k)
                    );

                //eta * x_h^2 * znt
                tmp4 = gsl_complex_mul_real(
                            gsl_matrix_complex_get(z,n,t)
                            ,gsl_matrix_get(eta,n,k)*gsl_matrix_get(x_h_abs2,k,t)
                        );

                tmp2.dat[0] += tmp3.dat[0] + tmp4.dat[0];
                tmp2.dat[1] += tmp3.dat[1] + tmp4.dat[1];
            }
            tmp2.dat[0] = tmp2.dat[0]/(double)N;
            tmp2.dat[1] = tmp2.dat[1]/(double)N;

            // printf("n = %d t = %d %f\n",n,t,tmp2.dat[0]);            
            //I_bから引く
            tmp3 = gsl_complex_sub(gsl_matrix_complex_get(I_b,n,t),tmp2);
            
            //I_bにset
            gsl_matrix_complex_set(I_b,n,t,tmp3);

        }
    }
    h_h_conju = GSLMatrixFree(h_h_conju);
}

//----------------------zの計算----------------------------------------
void culc_z()
{
    gsl_matrix_complex* tmp1;
    tmp1 = gsl_matrix_complex_calloc(N,T);
    gsl_matrix_complex_memcpy(tmp1,y);

    // y - I_b
    gsl_matrix_complex_sub(tmp1,I_b);

    //N0 + zeta
    gsl_matrix* tmp2;
    tmp2 = gsl_matrix_calloc(N,T);
    int n,t;

    for (n = 0; n < N; ++n)
    {
        for (t = 0; t < T; ++t)
        {
            double noise = N0;
            gsl_matrix_set(tmp2,n,t,noise);
        }
    }
    gsl_matrix_add(tmp2,zeta);

    //z = (y - I_b)/(N0 + zeta)
    for (n = 0; n < N; ++n)
    {
        for (t = 0; t < T; ++t)
        {
            gsl_complex tmp3 = gsl_complex_div_real(
                                gsl_matrix_complex_get(tmp1,n,t)
                                ,gsl_matrix_get(tmp2,n,t)
                                );
            gsl_matrix_complex_set(z,n,t,tmp3);
        }
    }
    tmp1 = GSLMatrixFree(tmp1);
    tmp2 = GSLRealMatrixFree(tmp2);

}