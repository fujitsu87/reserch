#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "../include/matrix.c"
#include "../include/random_number.c"
#include "../include/blas.c"
#include "../include/complex.c"
#include "Hadamard.c"
#include "../include/interleaver.c"
#include "init_functions.c"
#include "inspection.c"
#include <sys/stat.h>
#include <sys/types.h>

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
    free(user_p);
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
            if(tmp < 0)printf("%f\n",tmp );
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
            gsl_complex tmp3 = gsl_complex_div_real(
                                gsl_matrix_complex_get(tmp1,n,t)
                                            ,gsl_matrix_get(tmp2,n,t)
                                );
            gsl_matrix_complex_set(z,n,t,tmp3);
        }
    }
    tmp1 = GSLMatrixFree(tmp1);
    tmp2 = GSLRealMatrixFree(tmp2);

    // printf("%f\n",culc_abs2_all_element(N,T,z));
}
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

//----------------I_bの計算--------------------------------
void culc_I()
{

    int n,t,k;
    gsl_complex tmp;
    gsl_complex tmp2;
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
                //h_h^2 + xi * znt
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
                // if(tmp2 < 0){
                //     printf("exist - sec2 = %f\n",tmp2);
                //     tmp2=0.0;
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
            if(k < K/DK){
                tmp = user_p[k] - gsl_complex_abs2(gsl_matrix_complex_get(x_h,k,t));
            }
            else{
                tmp = user_p[k] - gsl_complex_abs2(gsl_matrix_complex_get(x_h,k,t));
            }

            // tmp = v_k(gsl_matrix_complex_get(x_h,k,t)
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
                tmp.dat[0] = (1.0 - a)*GSL_REAL(gsl_matrix_complex_get(x_h,k,t))
                            +a*fk(GSL_REAL(gsl_matrix_complex_get(x_b,k,t))
                                ,gsl_matrix_get(xi_b,k,t),k,t
                                );
                tmp.dat[1] = (1.0 - a)*GSL_IMAG(gsl_matrix_complex_get(x_h,k,t))
                            +a*fk(GSL_IMAG(gsl_matrix_complex_get(x_b,k,t))
                                ,gsl_matrix_get(xi_b,k,t),k,t
                                );
                gsl_matrix_complex_set(x_h,k,t,tmp);
            }

        }

    }

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

//------------------データ推定器----------------------------------
void data_estimation(FILE *fp_x)
{
    culc_zata();
    culc_I();
    culc_z();
   
    culc_xi_b();
    culc_x_b();

    culc_x_h();
    culc_xi();


    // output_estimatedata(fp_x);

}
//---------------------------------------------------------------
int main(int argc, char *argv[])
{
    if(argc != 3)
    {
        printf("Input error. Please input as below \n" );
        printf("./a.out SNRdB DIR\n");
        return 0;
    }

    FILE *fp_mse_h,*fp_mse_h_my,*fp_mse_h_other;
    FILE *fp_bit_err,*fp_bit_err_my,*fp_bit_err_other;
    FILE *fp_x;
    FILE *fp_sn_bit_err;

    char mse_h[100];
    char mse_h_my[100];
    char mse_h_other[100];
    char bit_err[100];
    char bit_err_my[100];
    char bit_err_other[100];
    char etc[100];
    char sn_bit_err[100];
    char output_path[100];
    
    int i,j,l,k;


    //mse
    //他セル　自セル含む
    double  *mse_h_n = (double *)malloc( sizeof(double) *BIG_LOOP*H_LOOP );
    double  *mse_h_en = (double *)malloc( sizeof(double) *ENSEMBLE);
    //自セル
    double  *mse_h_n_my = (double *)malloc( sizeof(double) *BIG_LOOP*H_LOOP );
    //他セル
    double  *mse_h_n_other = (double *)malloc( sizeof(double) *BIG_LOOP*H_LOOP );
    
    //BER
    //他セル　自セル含む
    double  *bit_err_n = (double *)malloc( sizeof(double) *BIG_LOOP*X_LOOP );
    double  *bit_err_en = (double *)malloc( sizeof(double) *ENSEMBLE );
    //自セル
    double  *bit_err_n_my = (double *)malloc( sizeof(double) *BIG_LOOP*X_LOOP );
    //他セル
    double  *bit_err_n_other = (double *)malloc( sizeof(double) *BIG_LOOP*X_LOOP );

    //I_b
    double  *I_b_n = (double *)malloc( sizeof(double) *BIG_LOOP*H_LOOP*X_LOOP );

    //zeta
    double  *zeta_n = (double *)malloc( sizeof(double) *BIG_LOOP*H_LOOP*X_LOOP );


    //乱数初期化
    RandomNumberInitialization(0);
   
    N0 = Pk/atof(argv[1]);
    sprintf(output_path,"./data/%s/",argv[2]);
    mkdir(output_path,0777);
    sprintf(mse_h,"./data/%s/mse_h_N%d_K%d_T%d_Tp%d_SN%.f.dat",argv[2],N,K,T,Tp,Pk/N0);
    sprintf(mse_h_my,"./data/%s/mse_h_my_N%d_K%d_T%d_Tp%d_SN%.f.dat",argv[2],N,K,T,Tp,Pk/N0);
    sprintf(mse_h_other,"./data/%s/mse_h_other_N%d_K%d_T%d_Tp%d_SN%.f.dat",argv[2],N,K,T,Tp,Pk/N0);
    sprintf(bit_err,"./data/%s/bit_err_N%d_K%d_T%d_Tp%d_SN%.f.dat",argv[2],N,K,T,Tp,Pk/N0);
    sprintf(bit_err_my,"./data/%s/bit_err_my_N%d_K%d_T%d_Tp%d_SN%.f.dat",argv[2],N,K,T,Tp,Pk/N0);    
    sprintf(bit_err_other,"./data/%s/bit_err_other_N%d_K%d_T%d_Tp%d_SN%.f.dat",argv[2],N,K,T,Tp,Pk/N0);
    sprintf(etc,"./data/%s/etc_N%d_K%d_T%d_Tp%d_SN%.f.dat",argv[2],N,K,T,Tp,Pk/N0);
    sprintf(sn_bit_err,"./data/%s/sn_bit_err_N%d_K%d_T%d_Tp%d.dat",argv[2],N,K,T,Tp);
    
    if ((fp_mse_h = fopen(mse_h, "w")) == NULL) {
        printf("can not open%s\n",mse_h);
        return 1;
    }
    if ((fp_mse_h_my = fopen(mse_h_my, "w")) == NULL) {
        printf("can not open%s\n",mse_h_my);
        return 1;
    }
    if ((fp_mse_h_other = fopen(mse_h_other, "w")) == NULL) {
        printf("can not open%s\n",mse_h_other);
        return 1;
    }
    if ((fp_bit_err = fopen(bit_err, "w")) == NULL) {
        printf("can not open %s\n",bit_err);
        return 1;
    }
    if ((fp_bit_err_my = fopen(bit_err_my, "w")) == NULL) {
        printf("can not open %s\n",bit_err_my);
        return 1;
    }
    if ((fp_bit_err_other = fopen(bit_err_other, "w")) == NULL) {
        printf("can not open %s\n",bit_err_other);
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

    clock_t start,end;
    start = clock();

    for (k = 0; k < ENSEMBLE; ++k)
    {
        init(atof(argv[1]));
        abs2_matrix(x_h,x_h_abs2,K,T);
        abs2_matrix(h_h,h_h_abs2,N,K);
        int count_h = 0;
        int count_x = 0;
        int count_all = 0;
        Repeat_flg = 0;

        for (i = 0; i < BIG_LOOP; ++i)
        {
            printf("i = %d H\n",i );

            //通信路推定
            init_eta()
            for (j = 0; j < H_LOOP; ++j)
            {
                if(j==0)gsl_matrix_set_zero(xi);
                // gsl_matrix_set_zero(eta);
                channel_estimation(fp_x);
                // fprintf(stdout, "%d %e\n",count_h,culc_complex_mse(N,K,h,h_h));
                mse_h_n[count_h] += culc_complex_mse(N,K,h,h_h)/(double)ENSEMBLE;
                mse_h_n_my[count_h] += culc_complex_mse_etc(N,K/DK,h,h_h,0)/(double)ENSEMBLE;
                mse_h_n_other[count_h] += culc_complex_mse_etc(N,K,h,h_h,K/DK)/(double)ENSEMBLE;
                ++count_h;
                abs2_matrix(h_h,h_h_abs2,N,K);

            
                // I_b_n[count_all] += culc_abs2_all_element(N,T,I_b)/(double)ENSEMBLE;
                // zeta_n[count_all] += culc_abs2_all_element(N,T,z)/(double)ENSEMBLE;
                ++count_all;
                // PrintMatrix(stdout,N,K,h_h);
                // PrintMatrix(stdout,N,K,h);

            }
           
            // PrintMatrix(fp_x,N,K,h_h);
            Repeat_flg = 1;

            gsl_matrix_complex_set_zero(z);
            gsl_matrix_complex_set_zero(I_b);
            gsl_matrix_set_zero(zeta);

            printf("i = %d X\n",i );
            //データ推定
            for (j = 0; j < X_LOOP; ++j)
            {
                data_estimation(fp_x);

                bit_err_n[count_x] += culc_bit_error_rate(fp_bit_err,0,K)/(double)ENSEMBLE;
                bit_err_n_my[count_x] += culc_bit_error_rate(fp_bit_err,0,K/DK)/(double)ENSEMBLE;
                bit_err_n_other[count_x] += culc_bit_error_rate(fp_bit_err,K/DK,K)/(double)ENSEMBLE;
                ++count_x;
                abs2_matrix(x_h,x_h_abs2,K,T);

                // I_b_n[count_all] += culc_abs2_all_element(N,T,I_b)/(double)ENSEMBLE;
                // zeta_n[count_all] += culc_abs2_all_element(N,T,z)/(double)ENSEMBLE;
                ++count_all;
                
            }
            
            gsl_matrix_complex_set_zero(z);
            gsl_matrix_complex_set_zero(I_b);
            gsl_matrix_set_zero(zeta);
        }
        mse_h_en[k] = culc_complex_mse(N,K,h,h_h);
        bit_err_en[k] = culc_bit_error_rate(fp_bit_err,0,K);

        finish();
    }

    end = clock();

    //ファイル書き込み
    for (i = 0; i < BIG_LOOP*H_LOOP; ++i)
    {
        fprintf(fp_mse_h,"%d %lf\n",i+1,mse_h_n[i]);
        fprintf(fp_mse_h_my,"%d %lf\n",i+1,mse_h_n_my[i]);
        fprintf(fp_mse_h_other,"%d %lf\n",i+1,mse_h_n_other[i]);
    }
    for (i = 0; i < BIG_LOOP*X_LOOP; ++i)
    {
        fprintf(fp_bit_err, "%d %lf\n",i+1,bit_err_n[i]);
        fprintf(fp_bit_err_my, "%d %lf\n",i+1,bit_err_n_my[i]);
        fprintf(fp_bit_err_other, "%d %lf\n",i+1,bit_err_n_other[i]);
    }

    //MSEとBER分散と標準偏差を求める
    double standard_deviation[2];
    culc_standard_deviation(standard_deviation,mse_h_n,mse_h_en,bit_err_n,bit_err_en);
    fprintf(fp_x,"MSE = %g \nmse_standard_deviation = %g\n",mse_h_n[BIG_LOOP*H_LOOP-1],standard_deviation[0]);
    fprintf(fp_x,"BER = %g \nber_standard_deviation = %g\n",bit_err_n[BIG_LOOP*X_LOOP-1],standard_deviation[1]);
    printf("BER = %g MSE = %g Computation time = %f[sec]\n",bit_err_n[BIG_LOOP*X_LOOP-1],mse_h_n[BIG_LOOP*H_LOOP-1],(double)(end-start)/ CLOCKS_PER_SEC);
    for(k=0;k<ENSEMBLE;k++)fprintf(fp_sn_bit_err,"%f %f\n",Pk/N0,bit_err_en[k]);

    create_plotfile(output_path);

    fclose(fp_mse_h);
    fclose(fp_mse_h_my);
    fclose(fp_mse_h_other);
    fclose(fp_bit_err);
    fclose(fp_x);
    fclose(fp_sn_bit_err);
    free(mse_h_n);
    free(bit_err_n); 
    free(mse_h_n_my);
    free(bit_err_n_my); 
    free(mse_h_n_other);
    free(bit_err_n_other); 
    free(mse_h_en);
    free(bit_err_en); 
    return 0;
}