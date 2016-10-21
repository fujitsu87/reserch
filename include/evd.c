#ifndef EVD_C
#define EVD_C 1

/* インクルードファイル */
#include <gsl/gsl_eigen.h>
#include "sort.c"
// 以下は上のファイル内でインクルード
//#include <stdio.h>
//#include <stdlib.h>
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_matrix.h>

/* プロトタイプ */
void EigenValue(gsl_matrix_complex *A, double lambda[]);
void EVD(gsl_matrix_complex *A, double lambda[], gsl_matrix_complex *U);

/* 関数本体 */
//----------------実固有値の計算----------------
//n次エルミート行列Aの固有値を計算し，n次配列lambda[]に格納する．
void EigenValue(gsl_matrix_complex *A, double lambda[])
{
	//行列保持用コピー
	gsl_matrix_complex *M = gsl_matrix_complex_alloc(A->size1, A->size2);
	gsl_matrix_complex_memcpy(M, A);

	//固有値計算
	{
		gsl_eigen_herm_workspace *w = gsl_eigen_herm_alloc(M->size1);
		gsl_vector *eval = gsl_vector_alloc(M->size1);

		gsl_eigen_herm(M, eval, w);

		//出力生成
		{
			int i;
			for(i = 0; i < eval->size; i++)
			{
				lambda[i] = gsl_vector_get(eval,i);
			}
			
			//ソート
			DSortDes(eval->size, lambda);
		}

		gsl_vector_free(eval);
		gsl_eigen_herm_free(w);
	}

	gsl_matrix_complex_free(M);
}

//----------------実固有値と固有ベクトルの計算----------------
//n次エルミート行列Aの固有値を計算し，n次配列lambda[]に格納する．また，固有ベクトルを計算し，n次正方行列Uに格納する．
void EVD(gsl_matrix_complex *A, double lambda[], gsl_matrix_complex *U)
{
	//行列保持用コピー
	gsl_matrix_complex *M = gsl_matrix_complex_alloc(A->size1, A->size2);
	gsl_matrix_complex_memcpy(M, A);

	//固有値計算
	{
		gsl_eigen_hermv_workspace *w = gsl_eigen_hermv_alloc(M->size1);
		gsl_vector *eval = gsl_vector_alloc(M->size1);

		gsl_eigen_hermv(M, eval, U, w);
		gsl_eigen_hermv_sort(eval, U, GSL_EIGEN_SORT_VAL_DESC);

		//出力生成
		{
			int i;
			for(i = 0; i < eval->size; i++)
			{
				lambda[i] = gsl_vector_get(eval, i);
			}
		}

		gsl_vector_free(eval);
		gsl_eigen_hermv_free(w);
	}

	gsl_matrix_complex_free(M);
}

#endif