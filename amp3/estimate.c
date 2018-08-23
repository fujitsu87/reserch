#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "../include/matrix.c"
#include "../include/vector.c"
#include "../include/random_number.c"
#include "../include/blas.c"
#include "../include/complex.c"
#include "Hadamard.c"
#include "../include/interleaver.c"
#include "init_functions.c"
#include "interference_functions.c"
#include "inspection.c"
#include "channel_functions.c"
#include "data_functions.c"
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
int Termination_d(gsl_matrix_complex* x_pre_one)
{
  int k,t;
  double eps;

  eps = 1.0e-2;
  for(k=0;k<K;k++) {
      for(t=0;t<T;t++) {
        gsl_complex t1 = gsl_matrix_complex_get(x_h,k,t);
        gsl_complex t2 = gsl_matrix_complex_get(x_h_pre,k,t);
        if(gsl_complex_abs(gsl_complex_sub(t1,t2))/gsl_complex_abs(t2)>eps) return 0;
      }
  }

  return 1;
}
//------------------終了処理-----------------------------------------
void finish()
{
    x = GSLMatrixFree(x); 
    y = GSLMatrixFree(y); 
    h = GSLMatrixFree(h); 
    w = GSLMatrixFree(w); 

    z_c = GSLMatrixFree(z_c); 
    z_d = GSLMatrixFree(z_d); 

    x_h = GSLMatrixFree(x_h);
    x_h_pre = GSLMatrixFree(x_h_pre);
    x_h_abs2 = GSLRealMatrixFree(x_h_abs2); 
    xi = GSLRealMatrixFree(xi);

    h_h = GSLMatrixFree(h_h);
    h_h_pre = GSLMatrixFree(h_h_pre);
    h_h_abs2 = GSLRealMatrixFree(h_h_abs2);
    eta = GSLRealMatrixFree(eta);

    zeta_c = GSLRealMatrixFree(zeta_c);
    zeta_d = GSLRealMatrixFree(zeta_d);

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
    
    int i,j,l,k,m;


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

    double **ber;
    ber = (double **)malloc(sizeof(double *) * BIG_LOOP);
    for (i=0;i<BIG_LOOP;i++) {
        ber[i] = (double *)malloc(sizeof(double) * X_LOOP);
    }
    for(i=0;i<BIG_LOOP;i++) {
        for(j=0;j<X_LOOP;j++) ber[i][j] = 0.0;
    }
    //自セル
    double  *bit_err_n_my = (double *)malloc( sizeof(double) *BIG_LOOP*X_LOOP );
    //他セル
    double  *bit_err_n_other = (double *)malloc( sizeof(double) *BIG_LOOP*X_LOOP );
    //zeta
    double  *zeta_n = (double *)malloc( sizeof(double) *BIG_LOOP*H_LOOP*X_LOOP );


    //乱数初期化
    RandomNumberInitialization(0);
   
    N0 = Pk/atof(argv[1]);
    sprintf(output_path,"./data/%s/",argv[2]);
    mkdir(output_path,0777);
    sprintf(mse_h,"./data/%s/mse_h_N%d_K%d_T%d_Tp%d_SN%.f.dat",argv[2],N,K,T,Tp,argv[1]);
    sprintf(mse_h_my,"./data/%s/mse_h_my_N%d_K%d_T%d_Tp%d_SN%.f.dat",argv[2],N,K,T,Tp,argv[1]);
    sprintf(mse_h_other,"./data/%s/mse_h_other_N%d_K%d_T%d_Tp%d_SN%.f.dat",argv[2],N,K,T,Tp,argv[1]);
    sprintf(bit_err,"./data/%s/bit_err_N%d_K%d_T%d_Tp%d_SN%.f.dat",argv[2],N,K,T,Tp,argv[1]);
    sprintf(bit_err_my,"./data/%s/bit_err_my_N%d_K%d_T%d_Tp%d_SN%.f.dat",argv[2],N,K,T,Tp,argv[1]);    
    sprintf(bit_err_other,"./data/%s/bit_err_other_N%d_K%d_T%d_Tp%d_SN%.f.dat",argv[2],N,K,T,Tp,argv[1]);
    sprintf(etc,"./data/%s/etc_N%d_K%d_T%d_Tp%d_SN%.f.dat",argv[2],N,K,T,Tp,argv[1]);
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

    gsl_matrix_complex* x_pre_one = gsl_matrix_complex_calloc(K,T);
    
    for (k = 0; k < ENSEMBLE; ++k)
    {
        init(atof(argv[1]));
        // init_xi_first();
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
            for (j = 0; j < H_LOOP; ++j)
            {

                channel_estimation(fp_x);

                mse_h_n[count_h] += culc_complex_mse(N,K,h,h_h)/(double)ENSEMBLE;
                mse_h_n_my[count_h] += culc_complex_mse_etc(N,K/DK,h,h_h,0)/(double)ENSEMBLE;
                mse_h_n_other[count_h] += culc_complex_mse_etc(N,K,h,h_h,K/DK)/(double)ENSEMBLE;
                ++count_h;
                abs2_matrix(h_h,h_h_abs2,N,K);
printf("MSE %d %g\n",i,culc_complex_mse(N,K,h,h_h));
                ++count_all;
            }
            gsl_matrix_complex_memcpy(h_h_pre,h_h); 
            Repeat_flg = 1;

            //データ推定
            gsl_matrix_complex_memcpy(x_pre_one,x_h);
            // gsl_matrix_complex_memcpy(h_h,h);
            // init_x_h();
            // init_xi_first();
            for (j = 0; j < T; ++j)
            {
                double ans[3];
                // data_estimation(fp_x);
                data_estimation(fp_x,j,ber[i]);
                // ans[0] = culc_bit_error_rate(fp_bit_err,0,K)/(double)ENSEMBLE;
                // ans[1] = culc_bit_error_rate(fp_bit_err,0,K/DK)/(double)ENSEMBLE;
                // ans[2] = culc_bit_error_rate(fp_bit_err,K/DK,K)/(double)ENSEMBLE;
                abs2_matrix(x_h,x_h_abs2,K,T);
              
//                 bit_err_n[count_x] += ans[0];
//                 bit_err_n_my[count_x] += ans[1];
//                 bit_err_n_other[count_x] += ans[2];
//                 ++count_x;
//                 ++count_all;
//                 gsl_matrix_complex_memcpy(x_pre_one,x_h);
            }
            for (m = 0; m < X_LOOP; ++m){
                printf("BER %d %g\n",i,ber[i][m]/(double)(k+1));
            }
            gsl_matrix_complex_memcpy(x_h_pre,x_h);
        }
        
        mse_h_en[k] = culc_complex_mse(N,K,h,h_h);
        bit_err_en[k] = culc_bit_error_rate(fp_bit_err,0,K);
        finish();
    }  
    x_pre_one =  GSLMatrixFree(x_pre_one);
    end = clock();

    //ファイル書き込み
    for (i = 0; i < BIG_LOOP*H_LOOP; ++i)
    {
        fprintf(fp_mse_h,"%d %lf\n",i+1,mse_h_n[i]);
        fprintf(fp_mse_h_my,"%d %lf\n",i+1,mse_h_n_my[i]);
        fprintf(fp_mse_h_other,"%d %lf\n",i+1,mse_h_n_other[i]);
    }
    // for (i = 0; i < BIG_LOOP*X_LOOP; ++i)
    // {
    //     fprintf(fp_bit_err, "%d %lf\n",i+1,bit_err_n[i]);
    //     fprintf(fp_bit_err_my, "%d %lf\n",i+1,bit_err_n_my[i]);
    //     fprintf(fp_bit_err_other, "%d %lf\n",i+1,bit_err_n_other[i]);
    // }
    for(i=0;i<BIG_LOOP;i++) {
        for(j=0;j<X_LOOP;j++) {
            ber[i][j] /= (double)ENSEMBLE;
            fprintf(fp_bit_err, "%d %lf\n",i*X_LOOP+j+1,ber[i][j]);
            // fprintf(fp_bit_err,"%g %g %d\n",dB,ber[i][j],i*iter_d+j+1);
        }
    }

    fprintf(fp_sn_bit_err,"%f %f\n",argv[1],ber[BIG_LOOP-1][X_LOOP-1]);
    //MSEとBER分散と標準偏差を求める
    double standard_deviation[2];
    culc_standard_deviation(standard_deviation,mse_h_n,mse_h_en,bit_err_n,bit_err_en);
    fprintf(fp_x,"MSE = %g \nmse_standard_deviation = %g\n",mse_h_n[BIG_LOOP*H_LOOP-1],standard_deviation[0]);
    fprintf(fp_x,"BER = %g \nber_standard_deviation = %g\n",bit_err_n[BIG_LOOP*X_LOOP-1],standard_deviation[1]);

    printf("BER = %g MSE = %g Computation time = %f[sec]\n",ber[BIG_LOOP-1][X_LOOP-1],mse_h_n[BIG_LOOP*X_LOOP-1],(double)(end-start)/ CLOCKS_PER_SEC);

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
    for (i = 0; i < BIG_LOOP; i++)
        free(ber[i]);
    free(ber);
    return 0;
}