#ifndef DCT_C
#define DCT_C 1

#include <gsl/gsl_vector.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <math.h>
#include "../include/fft.c"

/*----------------実ベクトルのDCT----------------*/
void DCT(gsl_vector *x, gsl_vector *y)
{
	int i;
    gsl_vector_complex *ext_x = gsl_vector_complex_alloc(2*x->size);    /*拡張した複素ベクトルの格納用*/
    gsl_vector_complex *ext_y = gsl_vector_complex_alloc(2*y->size);
    
    /*エラー処理*/
    if((x->size && y->size) != 1)
    {
        puts("void DCT(gsl_vector *x, gsl_vector *y): Vector length is invalid.");
        exit(0);
    }
    /*長さnの実ベクトルxを長さ2nの複素ベクトルext_xへ拡張*/
    for(i = 0;i < x->size;i++)
    {
        gsl_vector_complex_set(ext_x, i, gsl_complex_rect(gsl_vector_get(x, i), 0));
        gsl_vector_complex_set(ext_x, ext_x->size -1 - i, gsl_complex_rect(gsl_vector_get(x, i), 0));
    }
    /*離散フーリエ変換*/
    DFT(ext_x,ext_y);

    /*係数の処理および実ベクトルyへの復元*/
    gsl_vector_set(y, 0, GSL_REAL(gsl_complex_div_real(gsl_vector_complex_get(ext_y, 0), M_SQRT2)));
    for(i = 1;i < y->size;i++)
    {
         gsl_vector_set(y, i, GSL_REAL(gsl_complex_mul(gsl_vector_complex_get(ext_y, i), gsl_complex_polar(1, -i*M_PI_2/x->size))));
    }
    ext_x = GSLVectorFree(ext_x);
    ext_y = GSLVectorFree(ext_y);
}

/*----------------実ベクトルのIDCT----------------*/
void IDCT(gsl_vector *x, gsl_vector *y)
{
    int i;
    gsl_vector_complex *ext_x = gsl_vector_complex_alloc(2*x->size);    /*拡張した複素ベクトルの格納用*/
    gsl_vector_complex *ext_y = gsl_vector_complex_alloc(2*x->size);

    /*エラー処理*/
    if((x->size && y->size) != 1)
    {
        puts("void IDCT(gsl_vector *x, gsl_vector *y): Vector length is invalid.");
        exit(0);
    }
    /*係数の処理および複素ベクトルへの拡張*/
    gsl_vector_complex_set(ext_x, 0, gsl_complex_rect(gsl_vector_get(x, 0)*M_SQRT2, 0));
    for(i = 1;i < x->size;i++)
    {
         gsl_vector_complex_set(ext_x, i, gsl_complex_mul_real(gsl_complex_polar(1, i*M_PI_2/x->size), gsl_vector_get(x, i)));  
    }
    gsl_vector_complex_set(ext_x, x->size, gsl_complex_rect(0,0));
    for(i = x->size + 1;i < ext_x->size;i++)
    {
         gsl_vector_complex_set(ext_x, i, gsl_complex_conjugate(gsl_vector_complex_get(ext_x, ext_x->size - i)));
    }

    /*離散逆フーリエ変換*/
    IDFT(ext_x, ext_y);

    /*長さ2nの複素ベクトルext_yから長さnの実ベクトルyへ復元*/
    for(i = 0;i < y->size;i++)
    {
        gsl_vector_set(y, i, GSL_REAL(gsl_vector_complex_get(ext_y, i)));
    }
    ext_x = GSLVectorFree(ext_x);
    ext_y = GSLVectorFree(ext_y);
}
#endif
