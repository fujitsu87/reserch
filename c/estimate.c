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
#include "interference_functions.c"
#include "channel_functions.c"
#include "data_functions.c"
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
            printf("i = %d\n",i );

            //通信路推定
            gsl_matrix_complex_set_zero(h_h);
            init_eta();
            if(i==0)gsl_matrix_set_zero(xi);
            // if(i==0)init_xi_first();

            for (j = 0; j < H_LOOP; ++j)
            {

                channel_estimation(fp_x);

                mse_h_n[count_h] += culc_complex_mse(N,K,h,h_h)/(double)ENSEMBLE;
                mse_h_n_my[count_h] += culc_complex_mse_etc(N,K/DK,h,h_h,0)/(double)ENSEMBLE;
                mse_h_n_other[count_h] += culc_complex_mse_etc(N,K,h,h_h,K/DK)/(double)ENSEMBLE;
                ++count_h;
                abs2_matrix(h_h,h_h_abs2,N,K);

                printf("%f\n",culc_complex_mse(N,K,h,h_h));
            
                // I_b_n[count_all] += culc_abs2_all_element(N,T,I_b)/(double)ENSEMBLE;
                // zeta_n[count_all] += culc_abs2_all_element(N,T,z)/(double)ENSEMBLE;
                ++count_all;

            }


            Repeat_flg = 1;

            gsl_matrix_complex_set_zero(z);
            gsl_matrix_complex_set_zero(I_b);
            gsl_matrix_set_zero(zeta);

            //データ推定
            init_x_h();
            // gsl_matrix_set_zero(xi);
            init_xi_first();
            // PrintRealMatrix(stdout,N,K,eta);

            if(i==0)gsl_matrix_set_zero(eta);
            // if(i==0)tmp_set_eta(N,K,1.0);

            for (j = 0; j < X_LOOP; ++j)
            {
                data_estimation(fp_x);
                // fprintf(stdout, "%d %g\n",count_x,culc_bit_error_rate(fp_bit_err,0,K));
                bit_err_n[count_x] += culc_bit_error_rate(fp_bit_err,0,K)/(double)ENSEMBLE;
                bit_err_n_my[count_x] += culc_bit_error_rate(fp_bit_err,0,K/DK)/(double)ENSEMBLE;
                bit_err_n_other[count_x] += culc_bit_error_rate(fp_bit_err,K/DK,K)/(double)ENSEMBLE;
                ++count_x;
                abs2_matrix(x_h,x_h_abs2,K,T);

                // I_b_n[count_all] += culc_abs2_all_element(N,T,I_b)/(double)ENSEMBLE;
                // zeta_n[count_all] += culc_abs2_all_element(N,T,z)/(double)ENSEMBLE;
                ++count_all;
                
                // fprintf(stdout,"x_h = %g  xi_err = %g h_h = %g eta_err = %g\n",culc_complex_mse(K,T,x,x_h),fabs(culc_avr(K,T,xi)-culc_complex_mse(K,T,x,x_h)),culc_complex_mse(N,K,h,h_h),fabs(culc_avr(N,K,eta)-culc_complex_mse(N,K,h,h_h)));
            }
            
            gsl_matrix_complex_set_zero(z);
            gsl_matrix_complex_set_zero(I_b);
            gsl_matrix_set_zero(zeta);
            fprintf(fp_sn_bit_err,"%f %f\n",Pk/N0,culc_bit_error_rate(fp_bit_err,0,K));
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
    // for(k=0;k<ENSEMBLE;k++)fprintf(fp_sn_bit_err,"%f %f\n",Pk/N0,bit_err_en[k]);

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