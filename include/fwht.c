#ifndef FWHT_C
#define FWHT_C 1

#include "math.h"
#include "stdlib.h"
#include <gsl/gsl_vector.h>
#include "../include/vector.c"

/*プロトタイプ宣言*/
void FWHT(int m, gsl_vector *x, gsl_vector *y);
void IFWHT(int m, gsl_vector *x, gsl_vector *y);
void NonOrthogonalFWHT(int m, gsl_vector *x, gsl_vector *y);
void NonOrthogonalIFWHT(int m, gsl_vector *x, gsl_vector *y);

/*長さ2^mのベクトルxの直交FWHTを求め、結果をyに格納する関数*/
void FWHT(int m, gsl_vector *x, gsl_vector *y)
{
    int i;
    int n = (int)pow(2, m);

    NonOrthogonalFWHT(m, x, y);
    for(i = 0;i < n;i++)
        gsl_vector_set(y, i, gsl_vector_get(y,i)/sqrt(n));   
}
/*長さ2^mのベクトルxの直交IFWHTを求め、結果をyに格納する関数*/
void IFWHT(int m, gsl_vector *x, gsl_vector *y)
{
    FWHT(m, x, y);
}
/*長さ2^mのベクトルxの非直交FWHTを求め、結果をyに格納する関数*/
void NonOrthogonalFWHT(int m, gsl_vector *x, gsl_vector *y)
{  
    int i, n = (int)pow(2, m);

    
    if(m == 1)   /*ベクトルxの長さが2^1であれば*/
    {
        gsl_vector_set(y, 0, gsl_vector_get(x, 0) + gsl_vector_get(x, 1));
        gsl_vector_set(y, 1, gsl_vector_get(x, 0) - gsl_vector_get(x, 1));
    } 
    else
    {
        gsl_vector *hxu = gsl_vector_alloc(n/2); 
        gsl_vector *hxb = gsl_vector_alloc(n/2);
        gsl_vector_view xu, xb, yu, yb;
        xu = gsl_vector_subvector(x, 0, n/2);   /*ベクトルxの上半分をベクトルの像xuにとる*/
        xb = gsl_vector_subvector(x, n/2, n/2); /*ベクトルxの下半分をベクトルの像xbにとる*/
        yu = gsl_vector_subvector(y, 0, n/2);   /*ベクトルyの上半分をベクトルの像yuにとる*/
        yb = gsl_vector_subvector(y, n/2, n/2); /*ベクトルyの下半分をベクトルの像ybにとる*/
        NonOrthogonalFWHT(m-1, &xu.vector, hxu);/*ベクトルの像xuの非直交FWHTの結果をhxuへ格納*/
        NonOrthogonalFWHT(m-1, &xb.vector, hxb);/*ベクトルの像xbの非直交FWHTの結果をhxbへ格納*/
        gsl_vector_memcpy(&yu.vector, hxu); 
        gsl_vector_memcpy(&yb.vector, hxu);
        gsl_vector_add(&yu.vector, hxb);        /*ベクトル演算 yu = hxu + hxb*/
        gsl_vector_sub(&yb.vector, hxb);        /*ベクトル演算 yb = hxu - hxb*/
        hxu = GSLRealVectorFree(hxu);
        hxb = GSLRealVectorFree(hxb);
    }
}
/*長さ2^mのベクトルxの非直交IFWHTを求め、結果をyに格納する関数*/
void NonOrthogonalIFWHT(int m, gsl_vector *x, gsl_vector *y)
{  
    int i, n = (int)pow(2, m); 

    NonOrthogonalFWHT(m, x, y);
    for(i = 0;i < n;i++)
        gsl_vector_set(y, i, gsl_vector_get(y,i)/n);
}
#endif
