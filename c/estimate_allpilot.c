#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "../include/matrix.c"
#include "../include/random_number.c"
#include "../include/blas.c"
#include "init_functions_allpilot.c"
#include "inspection.c"


//----------------絶対値の2乗の行列を計算-------------------------------
void abs2_matrix(gsl_matrix_complex *x, gsl_matrix *y,int s1,int s2)
{
    int i,j;
    for (i = 0; i < s1; ++i)
    {
        for (j = 0; j < s2; ++j)
        {
            double tmp = gsl_complex_abs2(gsl_matrix_complex_get(x,i,j));
            gsl_matrix_set(y,i,j,tmp);

        }
    }
}
//------------------終了処理-----------------------------------------
void finish()
{
    x = GSLMatrixFree(x); 
    y = GSLMatrixFree(y); 
    h = GSLMatrixFree(h); 
    w = GSLMatrixFree(w); 

    z = GSLMatrixFree(z); 

    x_h = GSLMatrixFree(x_h);
    x_h_abs2 = GSLRealMatrixFree(x_h_abs2); 
    xi = GSLRealMatrixFree(xi);
    x_b = GSLMatrixFree(x_b);
    xi_b = GSLRealMatrixFree(xi_b);
    h_h = GSLMatrixFree(h_h);
    h_h_abs2 = GSLRealMatrixFree(h_h_abs2);
    eta = GSLRealMatrixFree(eta);
    I_b = GSLMatrixFree(I_b);
    zeta = GSLRealMatrixFree(zeta);

    pilot =  GSLRealMatrixFree(pilot);
}

//-----------zetaの計算-----------------------------------

void culc_zata()
{
    int n,t,k;

    abs2_matrix(x_h,x_h_abs2,K,T);
    abs2_matrix(h_h,h_h_abs2,N,K);

    double sec1,sec2,sec3;
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

            gsl_matrix_set(zeta,n,t,tmp);
        }
    }

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
            gsl_complex tmp3 = gsl_complex_div_real(gsl_matrix_complex_get(tmp1,n,t)
                                            ,gsl_matrix_get(tmp2,n,t));
            gsl_matrix_complex_set(z,n,t,tmp3);
        }
    }

    // PrintMatrix(stdout,N,T,z);
    tmp1 = GSLMatrixFree(tmp1);
    tmp2 = GSLRealMatrixFree(tmp2);
}
//---------------fk関数---------------------------------------------
double fk(double u,double v)
{
    double ans;
    double e1,e2,e3,t1,t2;
    e1 = exp(2*a_k*u/v);
    e2 = exp(-2*a_k*u/v);
    e3 = exp(a_k*a_k/v);

    t1 = a_k * (e1 - e2);
    t2 = e1 + e2 + 2 * (1/rho_k - 1) * e3;

    ans = t1/t2;


    return ans;
}

//------------------A関数-------------------------------------------
double diff_fk(double u,double v)
{
    double ans;
    double e1,e2,e3,t1,t2;
    e1 = exp(2*a_k*u/v);
    e2 = exp(-2*a_k*u/v);
    e3 = exp(a_k*a_k/v);

    t1 = 2*a_k/v*(e1*(a_k-1)+e2*(a_k-1));
    t2 = e1 + e2 + 2 * (1/rho_k - 1) * e3;

    t2 = t2 * t2;

    ans = t1/t2; 
    // printf("%d %f %f %f %f %e\n",Repeat_flg,u,v,t1,t2,ans);
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
                        ,v
                    );
        ans.dat[1] = GSL_IMAG(arg) 
                    * diff_fk(GSL_IMAG(gsl_matrix_complex_get(x_b,k,t))
                        ,v
                    );
        return ans;
    }
}

//----------------I_bの計算--------------------------------
void culc_I()
{
    int n,t,k;
    gsl_complex tmp;
    gsl_complex tmp2;
    gsl_complex tmp3;
    gsl_complex tmp4;
    abs2_matrix(x_h,x_h_abs2,K,T);

    gsl_matrix_complex *h_h_conju;
    h_h_conju = gsl_matrix_complex_calloc(K,N);
    ConjugateTranspose(h_h,h_h_conju);

    // PrintMatrix(stdout,K,N,h_h_conju);

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
                GSL_SET_COMPLEX(&tmp3,0,0);
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


            //I_bから引く
            tmp3 = gsl_complex_sub(gsl_matrix_complex_get(I_b,n,t),tmp2);
            //I_bにset
            gsl_matrix_complex_set(I_b,n,t,tmp3);

        }
    }

    h_h_conju = GSLMatrixFree(h_h_conju);
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
    for (n = 0; n < N; ++n)
    {
        for (k = 0; k < K; ++k)
        {
            GSL_SET_COMPLEX(&sec1,0,0);
            GSL_SET_COMPLEX(&sec2,0,0);
            GSL_SET_COMPLEX(&sec3,0,0);

            for (t = 0; t < T; ++t)
            {
                GSL_SET_COMPLEX(&tmp1,0,0);
                
                //h_h * x_h_conju　第1項の計算
                
                tmp1 = gsl_complex_mul(
                            gsl_matrix_complex_get(x_h_conju,t,k)
                            ,gsl_matrix_complex_get(z,n,t)
                        );
                sec1.dat[0] += tmp1.dat[0];
                sec1.dat[1] += tmp1.dat[1];
                
                //第3項の計算
                GSL_SET_COMPLEX(&tmp2,0,0);
                //Aktの計算
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

            }
            sec1.dat[0] = sec1.dat[0]/n_sqrt;
            sec1.dat[1] = sec1.dat[1]/n_sqrt;

            sec1 = gsl_complex_mul_real(sec1,gsl_matrix_get(eta,n,k));

            sec2 = gsl_complex_mul_real(
                    gsl_matrix_complex_get(h_h,n,k)
                    ,(1-gsl_matrix_get(eta,n,k))
                );


            sec3.dat[0] = sec3.dat[0]/(double)N;
            sec3.dat[1] = sec3.dat[1]/(double)N;
            sec3 = gsl_complex_mul_real(sec3,-gsl_matrix_get(eta,n,k));

            GSL_SET_COMPLEX(&tmp3,0,0);
            tmp3 = gsl_complex_add(sec1,sec2);
            tmp3 = gsl_complex_add(tmp3,sec3);
            gsl_matrix_complex_set(h_h,n,k,tmp3);
        }
    }

    x_h_conju = GSLMatrixFree(x_h_conju);
    h_h_conju = GSLMatrixFree(h_h_conju);


}
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
                tmp2 += gsl_complex_abs2(gsl_matrix_complex_get(h_h,n,k))/(N0 + gsl_matrix_get(zeta,n,t));
            }
            tmp1 = tmp2/(double)N;
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

                gsl_complex tmp1;
                tmp1.dat[0] = 0.0;
                tmp1.dat[1] = 0.0;

                double tmp2 = 0;
                gsl_complex tmp3,tmp4;

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
                // printf("%f\n",coef2*tmp2);
                tmp2 = 1 - coef2*tmp2;
                
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
        for (t = 0; t < T; ++t)
        {
            double tmp;
            tmp = Pk - gsl_complex_abs2(gsl_matrix_complex_get(x_h,k,t));
            // printf("%f %f\n",gsl_matrix_complex_get(x_h,k,t).dat[0],Pk);
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
        for (t = 0; t < T; ++t)
        {
            if(gsl_matrix_get(pilot,k,t) == 0)
            {
                tmp.dat[0] = fk(GSL_REAL(gsl_matrix_complex_get(x_b,k,t))
                                ,gsl_matrix_get(xi_b,k,t)
                                );
                tmp.dat[1] = fk(GSL_IMAG(gsl_matrix_complex_get(x_b,k,t))
                                ,gsl_matrix_get(xi_b,k,t)
                                );
                gsl_matrix_complex_set(x_h,k,t,tmp);
            }

        }

    }
}
//------------------通信路推定器----------------------------------
void channel_estimation()
{
    culc_zata();
    culc_I();
    culc_z();
    
    
    culc_eta();
    culc_h_h();
    
    /*
    gsl_matrix_set_zero(eta);
    gsl_matrix_complex_memcpy(h_h,h);
    */
}

//------------------データ推定器----------------------------------
void data_estimation()
{
    culc_zata();
    culc_I();
    culc_z();
    
    culc_xi_b();
    PrintMatrix(stdout,K,T,xi_b);
    culc_x_b();
    culc_xi();
    culc_x_h();
    // PrintRealMatrix(stdout,K,T,xi_b);
    // PrintMatrix(stdout,K,T,x_b);
}
//---------------------------------------------------------------
int main(int argc, char *argv[])
{
    FILE *fp_mse_h;
    FILE *fp_bit_err;
    FILE *fp_x;
    FILE *fp_sn_bit_err;

    char mse_h[100];
    char bit_err[100];
    char etc[100];
    char sn_bit_err[100];
    
    int i,j,l;
    int count_h = 1;
    int count_x = 1;
    init(atof(argv[1]));

    sprintf(mse_h,"./data_all_pilot/20161105_2/mse_h_N%d_K%d_T%d_Tp%d_SN%.f.dat",N,K,T,Tp,Pk/N0);
    sprintf(bit_err,"./data_all_pilot/20161105_2/bit_err_N%d_K%d_T%d_Tp%d_SN%.f.dat",N,K,T,Tp,Pk/N0);
    sprintf(etc,"./data_all_pilot/20161105_2/etc_N%d_K%d_T%d_Tp%d_SN%.f.dat",N,K,T,Tp,Pk/N0);
    sprintf(sn_bit_err,"./data_all_pilot/20161105_2/sn_bit_err_N%d_K%d_T%d_Tp%d.dat",N,K,T,Tp);
    
    if ((fp_mse_h = fopen(mse_h, "w")) == NULL) {
        printf("can not open%s\n",mse_h);
        return 1;
    }
    if ((fp_bit_err = fopen(bit_err, "w")) == NULL) {
        printf("can not open %s\n",bit_err);
        return 1;
    }
    if ((fp_x = fopen(etc, "w")) == NULL) {
        printf("can not open %s\n",etc);
        return 1;
    }
    if ((fp_sn_bit_err = fopen(sn_bit_err, "a")) == NULL) {
        printf("can not open %s\n",sn_bit_err);
        return 1;
    }

    // PrintMatrix(stdout,N,T,x_h);
    //x^2の計算
    abs2_matrix(x_h,x_h_abs2,K,T);
    abs2_matrix(h_h,h_h_abs2,N,K);

   

    clock_t start,end;
    start = clock();
    
    for (i = 0; i < 3; ++i)
    {
        
        //通信路推定
        for (j = 0; j < 100; ++j)
        {
            
            channel_estimation();
            fprintf(fp_mse_h, "%d %e\n",count_h,culc_complex_mse(N,K,h,h_h));
            ++count_h;
            // PrintMatrix(stdout,N,K,h_h)
            // PrintMatrix(stdout,N,T,I_b);;
            abs2_matrix(h_h,h_h_abs2,N,K);

        }
        PrintMatrix(stdout,N,K,h_h);
        // PrintMatrix(stdout,N,K,h);
        Repeat_flg = 1;


        //データ推定
        for (j = 0; j < 100; ++j)
        {


            data_estimation();
            fprintf(fp_bit_err, "%d %lf\n",count_x,culc_bit_error_rate());
            ++count_x;
            abs2_matrix(x_h,x_h_abs2,K,T);
            // PrintMatrix(stdout,K,T,x_h);
        }
        PrintMatrix(stdout,K,T,x_h);
        PrintMatrix(fp_x,K,T,x_h);
        PrintMatrix(fp_x,K,T,x);
    }
    end = clock();
    PrintRealMatrix(fp_x,K,T,pilot);
    // PrintRealMatrix(fp_x,K,T,xi);
    // PrintMatrix(stdout,N,K,h_h);
    // PrintMatrix(stdout,N,K,h);
    printf("BER = %g Computation time = %f[sec]\n",culc_bit_error_rate(),(double)(end-start)/ CLOCKS_PER_SEC);
    fprintf(fp_sn_bit_err,"%f %f\n",10*log10(Pk/N0),culc_bit_error_rate());

    finish();
    fclose(fp_mse_h);
    fclose(fp_bit_err);
    fclose(fp_x);
    fclose(fp_sn_bit_err);
    return 0;
}