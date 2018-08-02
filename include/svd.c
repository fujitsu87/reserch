#ifndef SVD_C
#define SVD_C 1

/* インクルードファイル */
#include "evd.c"
#include "blas.c"
#include <math.h>
// 以下は上のファイル内でインクルード
//#include <stdio.h>
//#include <stdlib.h>
//#include <gsl/gsl_complex.h>
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_matrix.h>


/* プロトタイプ */
void SingularValue(gsl_matrix_complex *A, double sigma[]);
unsigned int SVDRank(gsl_matrix_complex *A);
void SVDLeft(gsl_matrix_complex *A, double sigma[], gsl_matrix_complex *U);
void SVDRight(gsl_matrix_complex *A, double sigma[], gsl_matrix_complex *V);
void SVD(gsl_matrix_complex *A, double sigma[], gsl_matrix_complex *U, gsl_matrix_complex *V);

/* 関数本体 */
//---------------- 複素行列の非ゼロ特異値の計算 ----------------
// A : m×n行列，計算対象．
// sigma[] : k次配列，特異値格納用，ただしk=min(m,n)
void SingularValue(gsl_matrix_complex *A, double sigma[])
{
    //AAH, AHA の小さい方を選択
    gsl_matrix_complex *M;
    if(A->size1 < A->size2)
    {
        M = gsl_matrix_complex_alloc(A->size1, A->size1);
        AAH(A, M);
    }else{
        M = gsl_matrix_complex_alloc(A->size2, A->size2);
        AHA(A, M);
    }

    //固有値計算
    EigenValue(M, sigma);

    //出力生成
    {
        int i;
        for(i = 0; i < M->size1; i++)
        {
            double lambda = sigma[i];
            sigma[i] = (lambda < 0) ? 0 : sqrt(lambda);	//平方根代入
        }
    }

    gsl_matrix_complex_free(M);
}

//---------------- 複素行列のランクの計算 ----------------
// A : m×n行列，計算対象．
unsigned int SVDRank(gsl_matrix_complex *A)
{
    unsigned int k = (A->size1 < A->size2) ? A->size1 : A->size2;
    unsigned int rank;
    double *sigma = malloc(sizeof(double) * k);

    if(sigma == NULL){
        fprintf(stderr, "Error in '%s':", __func__);
        fprintf(stderr, "memory allocation error\n");
        exit(1);
    }

    //特異値計算
    SingularValue(A, sigma);

    //ゼロでない特異値の数を計測
    for(rank = 0; rank < k; rank++)
    {
        if(1E-16 >= sigma[rank])
            break;
    }

    //sigmaの解放
    free(sigma);

    return rank;
}

//---------------- 複素行列の特異値と左特異ベクトルの計算 ----------------
// A : m×n行列，計算対象．
// sigma[] : m次配列，特異値格納用．
// U : m×m複素正方行列，左特異ベクトル格納用．
void SVDLeft(gsl_matrix_complex *A, double sigma[], gsl_matrix_complex *U)
{
    //行列保持用コピー
    gsl_matrix_complex *_AAH_ = gsl_matrix_complex_alloc(A->size1, A->size1);
    AAH(A, _AAH_);

    //固有値分解
    EVD(_AAH_, sigma, U);

    //出力生成
    {
        int i;
        for(i = 0; i < _AAH_->size1; i++)
        {
            double lambda = sigma[i];
            if(lambda < 0)
            {
                sigma[i] = 0;								//特異値0
            }else{
                sigma[i] = sqrt(lambda);					//平方根代入
            }

        }
    }

    gsl_matrix_complex_free(_AAH_);
}


//---------------- 複素行列の特異値と右特異ベクトルの計算 ----------------
// A : m×n行列，計算対象．
// sigma[] : n次配列，特異値格納用．
// V : n×n複素正方行列，右特異ベクトル格納用．
void SVDRight(gsl_matrix_complex *A, double sigma[], gsl_matrix_complex *V)
{
    //行列保持用コピー
    gsl_matrix_complex *_AHA_ = gsl_matrix_complex_alloc(A->size2, A->size2);
    AHA(A, _AHA_);

    //固有値分解
    EVD(_AHA_, sigma, V);

    //出力生成
    {
        int i;
        for(i = 0; i < _AHA_->size1; i++)
        {
            double lambda = sigma[i];
            if(lambda < 0)
            {
                sigma[i] = 0;								//特異値0
            }else{
                sigma[i] = sqrt(lambda);					//平方根代入
            }
        }
    }

    gsl_matrix_complex_free(_AHA_);
}

//---------------- 複素行列の特異値と特異ベクトルの計算 ----------------
// A : m×n行列，計算対象．
// sigma[] : k次配列，特異値格納用，ただしk=min(m,n)
// U : m×m複素正方行列，左特異ベクトル格納用．
// V : n×n複素正方行列，右特異ベクトル格納用．
void SVD(gsl_matrix_complex *A, double sigma[], gsl_matrix_complex *U, gsl_matrix_complex *V)
{
    if(A->size1 < A->size2){
        double *s = malloc(sizeof(double) * A->size2);
        if(s == NULL){
            fprintf(stderr, "ERROR in '%s':", __func__);
            fprintf(stderr, "memory allocation is failed\n");
            exit(1);
        }

        SVDRight(A, s, V);

        {
            gsl_matrix_complex *AV = gsl_matrix_complex_alloc(A->size1, A->size2);
            gsl_matrix_complex *SI = gsl_matrix_complex_calloc(A->size2, A->size1);

            int i;
            for(i = 0; i < A->size1; i++){
                gsl_complex z = {(1.0 / s[i]), 0.0};
                gsl_matrix_complex_set(SI, i, i, z);
                sigma[i] = s[i];
            }

            AB(A, V, AV);
            AB(AV, SI, U);

            gsl_matrix_complex_free(AV);
            gsl_matrix_complex_free(SI);
        }

        free(s);
    }
    else
    {
        double *s = malloc(sizeof(double) * A->size1);
        if(s == NULL){
            fprintf(stderr, "ERROR in '%s':", __func__);
            fprintf(stderr, "memory allocation is failed\n");
            exit(1);
        }

        SVDLeft(A, s, U);

        {
            gsl_matrix_complex *AHU = gsl_matrix_complex_alloc(A->size2, A->size1);
            gsl_matrix_complex *SIH = gsl_matrix_complex_calloc(A->size1, A->size2);

            int i;
            for(i = 0; i < A->size2; i++){
                gsl_complex z = {(1.0 / s[i]), 0.0};
                gsl_matrix_complex_set(SIH, i, i, z);
                sigma[i] = s[i];
            }

            AHB(A, U, AHU);
            AB(AHU, SIH, V);

            gsl_matrix_complex_free(AHU);
            gsl_matrix_complex_free(SIH);
        }

        free(s);
    }

}

#endif