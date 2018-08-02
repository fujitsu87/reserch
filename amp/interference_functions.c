void abs2_matrix(gsl_matrix_complex *x, gsl_matrix *y,int s1,int s2);

//-----------zetaの計算-----------------------------------
//flg = 0 channel flg = 1 data
void culc_zata(int flg)
{
    gsl_matrix *zeta;
    if (flg == 0) 
    {
        zeta = zeta_c;
    }
    else
    {
        zeta = zeta_d;
    }
    int n,t,k;

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
                sec2 = gsl_matrix_get(eta,n,k) * gsl_complex_abs2(gsl_matrix_complex_get(x_h,k,t));
                sec3 = gsl_complex_abs2(gsl_matrix_complex_get(h_h,n,k)) * gsl_matrix_get(xi,k,t);
                tmp += sec1 + sec2 + sec3;
            }
            tmp = tmp/(double)N + N0;
            if(tmp < 0)printf("%f\n",tmp);
            // else if(tmp > 1)tmp = 1.0;
            gsl_matrix_set(zeta,n,t,tmp);
        }
    }
}

//----------------------zの計算----------------------------------------
//flg = 0 channel flg = 1 data
void culc_z(int flg)
{
    gsl_matrix_complex *z;
    gsl_matrix *zeta;
    if (flg == 0) 
    {
        z = z_c;
        zeta = zeta_c;
    }
    else
    { 
        z = z_d;
        zeta = zeta_d;
    }
    int n,t,k;
    for( n = 0 ; n < N ; ++n )
    {
        for( t = 0 ; t < T ; ++t )
        {
            gsl_complex sec1,sec2,sec3,sec4,ans;
            
            // y
            sec1 = gsl_matrix_complex_get(y,n,t);
            // sec1 = gsl_complex_div_real(
            //     sec1,gsl_matrix_get(zeta,n,t)
            // );
            
            // h_h * x_h
            sec2.dat[0] = 0.0;
            sec2.dat[1] = 0.0;
            for( k = 0 ; k < K ; ++k )
            {
                sec2 = gsl_complex_add(
                    gsl_complex_mul(
                        gsl_matrix_complex_get(h_h,n,k),
                        gsl_matrix_complex_get(x_h,k,t)
                    ),
                    sec2
                );
            }
            sec2 = gsl_complex_mul_real(sec2,-1.0);
            sec2 = gsl_complex_div_real(sec2,sqrt((double)N)); 

            //(eta * x_h_pre * x_h * z^c)/zeta_c
            sec3.dat[0] = 0.0;
            sec3.dat[1] = 0.0;
            for( k = 0 ; k < K ; ++k )
            {
                //eta * x_h_pre
                gsl_complex tmp1 = gsl_complex_mul_real(
                    gsl_complex_conjugate(
                        gsl_matrix_complex_get(x_h_pre,k,t)
                    ),
                    gsl_matrix_get(eta,n,k)
                );
                //x_h * z_c
                gsl_complex tmp2 = gsl_complex_mul(
                    gsl_matrix_complex_get(x_h,k,t),
                    gsl_matrix_complex_get(z_c,n,t)
                );
                sec3 = gsl_complex_add(
                    gsl_complex_div_real(
                        gsl_complex_mul(
                            tmp1,
                            tmp2
                        ),
                        gsl_matrix_get(zeta_c,n,t)
                    ),
                    sec3
                );
            }
            sec3 = gsl_complex_div_real(sec3,(double)N); 

            //(h_h_pre * h_h * xi * z_d)/zeta_d
            sec4.dat[0] = 0.0;
            sec4.dat[1] = 0.0;
            for( k = 0 ; k < K ; ++k )
            {
                //h_h_pre * h_h
                gsl_complex tmp1 = gsl_complex_mul(
                    gsl_complex_conjugate(
                        gsl_matrix_complex_get(h_h_pre,n,k)
                    ),
                    gsl_matrix_complex_get(h_h,n,k)
                );
                //xi * z_d
                gsl_complex tmp2 = gsl_complex_mul_real(
                    gsl_matrix_complex_get(z_d,n,t),
                    gsl_matrix_get(xi,k,t)
                );
                sec4 = gsl_complex_add(
                    gsl_complex_div_real(
                        gsl_complex_mul(
                            tmp1,
                            tmp2
                        ),
                        gsl_matrix_get(zeta_d,n,t)
                    ),
                    sec4
                );
            }
            sec4 = gsl_complex_div_real(sec4,(double)N); 
            ans = gsl_complex_add(sec1,sec2);
            ans = gsl_complex_add(ans,sec3);
            ans = gsl_complex_add(ans,sec4);
            gsl_matrix_complex_set(z,n,t,ans);
        }
    }
}

