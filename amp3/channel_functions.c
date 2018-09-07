//------------------culc_eta_b------------------------------------
void culc_eta()
{
    int n,k,t;
    double tmp;

    for (n = 0; n < N; ++n)
    {
        for (k = 0; k < K; ++k)
        {
            tmp = 0.0;
            for (t = 0; t < T; ++t)
            {
                tmp += gsl_matrix_get(x_h_abs2,k,t)/(gsl_matrix_get(zeta_c,n,t));
            }
            // tmp = tmp/(double)N; 
            // tmp = 1.0/((double)N*(1.0 + tmp));
            tmp = 1.0/tmp; 
            tmp = tmp/(1.0 + (double)N*tmp);
            gsl_matrix_set(eta,n,k,tmp);
        }
    }
    // PrintRealMatrix(stdout,N,T,zeta_c);
}

void culc_h_h()
{
    int n,k,t;

    for (n = 0; n < N; ++n)
    {
        for (k = 0; k < K; ++k)
        {
            gsl_complex sec1,sec2,sec3,ans;

            // h_h
            sec1 = gsl_complex_mul_real(
                gsl_matrix_complex_get(h_h,n,k),
                1 - (double)N*gsl_matrix_get(eta,n,k)
            );

            //eta_b * sum(x_h*z/zeta_c)
            sec2.dat[0] = 0.0;
            sec2.dat[1] = 0.0;
            for (t = 0; t < T; ++t)
            {
                sec2 = gsl_complex_add(
                    gsl_complex_div_real(
                        gsl_complex_mul(
                            gsl_complex_conjugate(
                                gsl_matrix_complex_get(x_h,k,t)
                            ),
                            gsl_matrix_complex_get(z_c,n,t)
                        ),
                        gsl_matrix_get(zeta_c,n,t)
                    )
                    ,sec2
                );
            }
            sec2 = gsl_complex_mul_real(
                sec2,
                gsl_matrix_get(eta,n,k)
            );
            // sec2 = gsl_complex_div_real(sec2,sqrt((double)N)); 
            //- eta_b * h_h_pre * sum(xi*z_d*z/zeta_d*zeta)
            sec3.dat[0] = 0.0;
            sec3.dat[1] = 0.0;
            for (t = 0; t < T; ++t)
            {
                gsl_complex tmp;
                //xi * z_d_con
                tmp = gsl_complex_mul_real(
                    gsl_complex_conjugate(
                        gsl_matrix_complex_get(z_d,n,t)
                    ),
                    gsl_matrix_get(xi,k,t)
                );
                //xi * z_d_con * z_c
                tmp = gsl_complex_mul(
                    tmp,
                    gsl_matrix_complex_get(z_c,n,t)
                );
                tmp =  gsl_complex_div_real(
                    tmp,
                    gsl_matrix_get(zeta_d,n,t)*gsl_matrix_get(zeta_c,n,t)
                );
                sec3 = gsl_complex_add(
                    tmp,
                    sec3
                );
            }
            // sec3 = gsl_complex_div_real(sec3,(double)N); 
            sec3 = gsl_complex_mul(
                sec3,
                gsl_matrix_complex_get(h_h,n,k)
            );
            sec3 = gsl_complex_mul_real(
                sec3,
                gsl_matrix_get(eta,n,k)
            );
            sec3 = gsl_complex_mul_real(
                sec3,
                -1.0
            );

            ans = gsl_complex_add(
                sec1,sec2
            );
            ans = gsl_complex_add(
                ans,sec3
            );
            gsl_complex pre = gsl_matrix_complex_get(h_h,n,k);
            ans.dat[0] = a_h * ans.dat[0] + (1 - a_h) * pre.dat[0];
            ans.dat[1] = a_h * ans.dat[1] + (1 - a_h) * pre.dat[1];
            gsl_matrix_complex_set(h_h,n,k,ans);
        }
    }
}

//------------------通信路推定器----------------------------------
void channel_estimation(FILE *fp_x)
{
    int n,k,t;
    for(n=0;n<N;n++) {
        for(t=0;t<T;t++) {
            gsl_complex z = Computation_z(n,t,x_h,h_h);
            gsl_matrix_complex_set(z_c,n,t,z);
            double zeta = Computation_zeta(n,t);
            gsl_matrix_set(zeta_c,n,t,zeta);
        }
    }
    
    culc_eta();
    culc_h_h();
}
