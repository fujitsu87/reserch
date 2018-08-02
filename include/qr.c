#ifndef QRV_C
#define QRV_C 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>		//GSL_MINで使用
#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>


//プロトタイプ宣言
void QRV(unsigned int m, unsigned int n, gsl_matrix_complex *A, gsl_matrix_complex *Q, gsl_matrix_complex *R);

 
// Function: _QR_decomp_complex
// ハウスホルダ変換を行うサブルーチン
// gsl_linalg_QR_decomp関数を参考にして複素数拡張
int QR_decomp_complex(unsigned int m, unsigned int n, gsl_matrix_complex *a, gsl_vector_complex *tau){

	if (m < n){
		GSL_ERROR (" Matrix size error M >= N ", GSL_EBADLEN);

    }else if (tau->size != GSL_MIN (m, n)){
        GSL_ERROR ("R matrix must be M x N", GSL_ENOTSQR);

	}else{
		size_t i = 0;

        for (i = 0; i < GSL_MIN (m, n); i++) {
            // i列ベクトルの像
            gsl_vector_complex_view  c_full = gsl_matrix_complex_column (a, i);

            // 列ベクトルのi列以降の像を取得
            gsl_vector_complex_view  c = gsl_vector_complex_subvector (&(c_full.vector), i, m - i);

            //  ハウスホルダ変換
            gsl_complex tau_i = gsl_linalg_complex_householder_transform (&(c.vector));

            gsl_vector_complex_set (tau, i, tau_i);

            //i 行 i + 1 列以降を Householder変換
            if (i + 1 < n) {
                gsl_matrix_complex_view b = gsl_matrix_complex_submatrix (a, i, i + 1, m - i, n - (i + 1));
                gsl_linalg_complex_householder_hm (tau_i, &(c.vector), &(b.matrix));
            }
        }
    }
	return GSL_SUCCESS;
}


// Function: QR_unpack_complex
// ハウスホルダ変換後の行列からQとRを取り出すサブルーチン
// gsl_linalg_QR_unpack関数を参考にして複素数拡張
int QR_unpack_complex(const gsl_matrix_complex * QR, const gsl_vector_complex * tau, gsl_matrix_complex * Q, gsl_matrix_complex * R){

    const size_t M = QR->size1;
    const size_t N = QR->size2;

    if (Q->size1 != M || Q->size2 != M){
        GSL_ERROR ("Q matrix must be M x M", GSL_ENOTSQR);
    
    }else if (R->size1 != M || R->size2 != N){
        GSL_ERROR ("R matrix must be M x N", GSL_ENOTSQR);
    
    }else if (tau->size != GSL_MIN (M, N)){
        GSL_ERROR ("size of tau must be MIN(M,N)", GSL_EBADLEN);
    
    }else{
        size_t i, j;

        gsl_matrix_complex_set_identity (Q);

        for (i = GSL_MIN (M, N); i-- > 0;){
            gsl_vector_complex_const_view c = gsl_matrix_complex_const_column (QR, i);
            gsl_vector_complex_const_view h = gsl_vector_complex_const_subvector (&c.vector,i, M - i);
            gsl_matrix_complex_view m = gsl_matrix_complex_submatrix (Q, i, i, M - i, M - i);
            gsl_complex ti = gsl_vector_complex_get (tau, i);
            gsl_linalg_complex_householder_hm (ti, &h.vector, &m.matrix);
        }

        for (i = 0; i < M; i++){
            for (j = 0; j < i && j < N; j++) gsl_matrix_complex_set (R, i, j, gsl_complex_rect(0.0,0.0));
            for (j = i; j < N; j++) gsl_matrix_complex_set (R, i, j, gsl_matrix_complex_get (QR, i, j));
        }

        return GSL_SUCCESS;
    }
}


void QRV(unsigned int m, unsigned int n, gsl_matrix_complex *A, gsl_matrix_complex *Q, gsl_matrix_complex *R){

    gsl_vector_complex *tau = gsl_vector_complex_alloc (GSL_MIN (m, n));

    // ハウスホルダ変換
    QR_decomp_complex (m,n,A,tau);

    #ifdef DEBUG_QRV
        printf("Debug: Matrix QR --------------------\n\n");
        PrintMatrix(stdout,3,3,A);
    #endif

    // QとRを取り出し
    QR_unpack_complex(A, tau, Q, R);

}

#endif