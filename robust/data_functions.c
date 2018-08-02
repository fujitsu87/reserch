//---------------fk関数---------------------------------------------
double fk(double u,double v,int k,int t)
{
    double ans;
    double e1,e2,e3,t1,t2;
    // if(k < K/DK){
        e1 = exp(2*a_k*u/v);
        e2 = exp(-2*a_k*u/v);
        e3 = exp(a_k*a_k/v);

        t1 = a_k * (e1 - e2);
    //}
    // else{
    //     e1 = exp(2*sqrt(cell_diff)*a_k*u/v);
    //     e2 = exp(-2*sqrt(cell_diff)*a_k*u/v);
    //     e3 = exp(0.5*a_k*a_k/v);

    //     t1 = sqrt(cell_diff)*a_k * (e1 - e2);
    // }
    t2 = e1 + e2 + 2 * (1/rho_k - 1) * e3;

    ans = t1/t2;

    return ans;
}

//------------------v_k関数-------------------------------------------
double v_k(gsl_complex u,double v,int k,int t)
{
    double ans;
    double e1,e2,e3,t1,t2,fk_real,fk_imag;
    double u2;
    u2 = gsl_complex_abs2(u);

    e1 = exp(a_k*u2/v);
    e2 = exp((-1)*a_k*u2/v);
    e3 = exp(a_k*a_k/(2*v));

    t1 = a_k * (e1 - e2);
    t2 = e1 + e2 + 2 * (1/rho_k - 1) * e3;

    fk_real = fk(GSL_REAL(u),v,k,t);
    fk_imag = fk(GSL_IMAG(u),v,k,t);


    ans = t1/t2 - (fk_real*fk_real + fk_imag*fk_imag);

    return ans;
}

//------------------A関数-------------------------------------------
double diff_fk(double u,double v,int k,int t)
{
    double ans;
    double e1,e2,e3,t1,t2;
    // if(k < K/DK){
        e1 = exp(2*a_k*u/v);
        e2 = exp(-2*a_k*u/v);
        e3 = exp(a_k*a_k/v);
        t1 = 2 * a_k* a_k / v * (4 + 2 * (1/rho_k - 1)* e3 *(e1 + e2) );
    // }
    /*else{
        e1 = exp(2*sqrt(cell_diff)*a_k*u/v);
        e2 = exp(-2*sqrt(cell_diff)*a_k*u/v);
        e3 = exp(cell_diff*a_k*a_k/v);
        t1 = 2 *cell_diff* a_k* a_k / v * (4 + 2 * (1/rho_k - 1)* e3 *(e1 + e2) );
    }*/
    t2 = e1 + e2 + 2 * (1/rho_k - 1) * e3;

    t2 = t2 * t2;

    ans = t1/t2; 

    return ans;
}


gsl_complex A(gsl_complex arg,int k,int t)
{
    gsl_complex ans;

    if(Repeat_flg == 0)
    {
        ans.dat[0] = 0.0;
        ans.dat[1] = 0.0;
        return ans;
    }
    else {
        double v = gsl_matrix_get(xi_b,k,t);
        ans.dat[0] = GSL_REAL(arg) 
                    * diff_fk(GSL_REAL(gsl_matrix_complex_get(x_b,k,t))
                        ,v,k,t
                    );
        ans.dat[1] = GSL_IMAG(arg) 
                    * diff_fk(GSL_IMAG(gsl_matrix_complex_get(x_b,k,t))
                        ,v,k,t
                    );
        return ans;
    }
}

//------------------culc_xi_b------------------------------------
void culc_xi_b()
{
    int k,t,n;

    for (k = 0; k < K; ++k)
    {
        for (t = 0; t < T; ++t)
        {
            double tmp1 = 0;
            double tmp2 = 0;
            for (n = 0; n < N; ++n)
            {
                // tmp2 += gsl_complex_abs2(gsl_matrix_complex_get(h_h,n,k))/(N0 + gsl_matrix_get(zeta,n,t));
                // tmp2 += gsl_complex_abs2(gsl_matrix_complex_get(h_h,n,k))/(N0);
                tmp2 += gsl_complex_abs2(gsl_matrix_complex_get(h_h,n,k))/gsl_matrix_get(robust_matrix,n,t);
            }
            tmp1 = tmp2/(double)N;
            if(tmp1 < 0)
            {
                printf("%f\n",tmp1);
            }
            gsl_matrix_set(xi_b,k,t,1/tmp1);
        }
    }
}

//------------------culc_x_b------------------------------------
void culc_x_b()
{
    int k,t,n;
    double tmp1,tmp2;

    for (k = 0; k < K; ++k)
    {
        for (t = 0; t < T; ++t)
        {
            if(gsl_matrix_get(pilot,k,t) == 0)
            {
                //sec1 ... 第一項　sec2 ... 第2項
                gsl_complex sec1,sec2;

                //tmp1 ...  sum (h_h* * z )
                gsl_complex tmp1;
                tmp1.dat[0] = 0.0;
                tmp1.dat[1] = 0.0;

                //tmp2 ... 1 - coef2*sum (eta*z^2)
                double tmp2 = 0;
                
                // tmp3 ... h_h* tmp4...z
                gsl_complex tmp3,tmp4;
                
                //tmp5 ... eta  tmp6...|z|^2
                double tmp5,tmp6;

                double coef1 = gsl_matrix_get(xi_b,k,t)/sqrt((double)N);
                double coef2 = gsl_matrix_get(xi_b,k,t)/(double)N;
                for (n = 0; n < N; ++n)
                {
                    tmp3 = gsl_complex_conjugate(gsl_matrix_complex_get(h_h,n,k));
                    tmp4 = gsl_matrix_complex_get(z,n,t);
                    tmp1 = gsl_complex_add(gsl_complex_mul(tmp3,tmp4),tmp1);

                    tmp5 = gsl_matrix_get(eta,n,k);
                    tmp6 = gsl_complex_abs2(gsl_matrix_complex_get(z,n,t));
                    tmp2 += tmp5*tmp6;
                }
                sec1 = gsl_complex_mul_real(tmp1,coef1);


                // coef2 = 0;


                tmp2 = 1 - coef2*tmp2;
                // if(tmp2 < 0){
                    // printf("exist - sec2 = %f\n",tmp2);
                    // tmp2=0.0;
                // }
                sec2 = gsl_complex_mul_real(gsl_matrix_complex_get(x_h,k,t),tmp2);
                gsl_matrix_complex_set(x_b,k,t,gsl_complex_add(sec1,sec2));

            }
        }
    }

}

//------------------culc_xi------------------------------------
void culc_xi()
{
    int k,t;
    for (k = 0; k < K; ++k)
    {
        a_k = sqrt(user_p[k]/(2*rho_k));
        for (t = 0; t < T; ++t)
        {
            double tmp;
            //旧ver 2017.11.14
            tmp = user_p[k] - gsl_complex_abs2(gsl_matrix_complex_get(x_h,k,t));

            // tmp = v_k(gsl_matrix_complex_get(x_b,k,t)
            //         ,gsl_matrix_get(xi_b,k,t),k,t);

            if(tmp < 0.0)
            {
                tmp=0.0;
            }
            gsl_matrix_set(xi,k,t,tmp);
        }
    }
}
//------------------culc_x_h------------------------------------
void culc_x_h()
{
    int k,t;
    gsl_complex tmp;
    for (k = 0; k < K; ++k)
    {
        a_k = sqrt(user_p[k]/(2*rho_k));
        for (t = 0; t < T; ++t)
        {
            if(gsl_matrix_get(pilot,k,t) == 0)
            {
                tmp.dat[0] = (1.0 - a_x)*GSL_REAL(gsl_matrix_complex_get(x_h,k,t))
                            +a_x*fk(GSL_REAL(gsl_matrix_complex_get(x_b,k,t))
                                ,gsl_matrix_get(xi_b,k,t),k,t
                                );
                tmp.dat[1] = (1.0 - a_x)*GSL_IMAG(gsl_matrix_complex_get(x_h,k,t))
                            +a_x*fk(GSL_IMAG(gsl_matrix_complex_get(x_b,k,t))
                                ,gsl_matrix_get(xi_b,k,t),k,t
                                );
                gsl_matrix_complex_set(x_h,k,t,tmp);
            }

        }

    }

}


//------------------データ推定器----------------------------------
void data_estimation(FILE *fp_x)
{
    culc_robust();
    // culc_zata();
    culc_I();
    culc_z();
   
    culc_xi_b();
    culc_x_b();

    culc_x_h();
    culc_xi();

    // output_estimatedata(fp_x);

}