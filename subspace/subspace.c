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
#include "svd_subspace.c"
#include <sys/stat.h>
#include <sys/types.h>


//------------------終了処理-----------------------------------------
void finish()
{
    x = GSLMatrixFree(x); 
    y = GSLMatrixFree(y); 
    h = GSLMatrixFree(h); 
    w = GSLMatrixFree(w); 

    z = GSLMatrixFree(z); 

    x_h = GSLMatrixFree(x_h);

    h_h = GSLMatrixFree(h_h);

    pilot =  GSLRealMatrixFree(pilot);

    s =  GSLMatrixFree(s);
    y_sub =  GSLMatrixFree(y_sub);
    h_sub =  GSLMatrixFree(h_sub);
    h_sub =  GSLMatrixFree(h_sub_true);
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

    int i,j,l,k;
    double mse_result = 0;

    FILE *fp_sn_mes;
    char sn_mes[100];
    sprintf(sn_mes,"./data/%s/sn_mes_N%d_K%d_T%d_Tp%d.dat",argv[2],N,K,T,Tp);
    if ((fp_sn_mes = fopen(sn_mes, "a")) == NULL) {
        printf("can not open %s\n",sn_mes);
        return 1;
    }
    
    
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


    clock_t start,end;
    start = clock();

    for (k = 0; k < ENSEMBLE; ++k)
    {
        init(atof(argv[1]));

        int count_h = 0;
        int count_x = 0;
        int count_all = 0;
        Repeat_flg = 0;
        

        // PrintMatrix(stdout,K,T,x_h);
        //通信路推定
        //sを計算
        svd();

        //y_subを計算
        AHB(s,y,y_sub);
        // fprintf(stdout, "\n\n---y_sub---\n");
        // PrintMatrix(stdout, K, T, y_sub);

        //h_subを計算
        ABH(y_sub,x_h,h_sub);
        gsl_complex tmp1;
	    GSL_SET_COMPLEX(&tmp1,1/sqrt(N),0);
	    gsl_matrix_complex_scale(h_sub,tmp1);
        // fprintf(stdout, "\n\n---h_sub---\n");
        // PrintMatrix(stdout, K/2, K, h_sub);

        //h_subの真値を計算
	    AHB(s,h,h_sub_true);
        GSL_SET_COMPLEX(&tmp1,1/sqrt(N),0);
	    gsl_matrix_complex_scale(h_sub_true,tmp1);
        // fprintf(stdout, "\n\n---h_sub_true---\n");
        // PrintMatrix(stdout, K/2, K/2, h_sub_true);


        // printf("mse = %g\n",culc_complex_mse(K/2,K,h_sub,h_sub_true));
        //データ推定
        
        // //データ書き込み
        mse_result = mse_result + culc_complex_mse(K/2,K,h_sub,h_sub_true);

        finish();
    }
    fprintf(fp_sn_mes,"%lf %lf\n",Pk/N0,mse_result/ENSEMBLE);
    end = clock();

    fclose(fp_sn_mes);
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