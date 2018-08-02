#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include "../include/matrix.c"
#include "../include/random_number.c"
#include "../include/blas.c"
#include "../include/complex.c"
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>


//----------------N0+ \zeta_nt^{i}-------------------------------
void culc_data_setting1(double *data1)
{
    int n,t;
    int index = 0;
    for(n=0;n<N;++n)
    {
        for(t=0;t<Tp;++t)
        {
            data1[index]= gsl_complex_abs(
                gsl_complex_sub(
                    gsl_matrix_complex_get(y,n,t),
                    gsl_matrix_complex_get(I_b,n,t)
                )
            );
            ++index;
        }
    }
    gsl_sort(data1, 1, N*Tp);
}
void culc_data_setting2(double *data2)
{
    int n,t;
    int index = 0;
    for(n=0;n<N;++n)
    {
        for(t=0;t>=Tp && t<(T-Tp);++t)
        {
            data2[index]= gsl_complex_abs(
                gsl_complex_sub(
                    gsl_matrix_complex_get(y,n,t),
                    gsl_matrix_complex_get(I_b,n,t)
                )
            );
            ++index;
        }
    }
    gsl_sort(data2, 1, N*(T - 2*Tp));
}
void culc_data_setting3(double *data3)
{
    int n,t;
    int index = 0;
    for(n=0;n<N;++n)
    {
        for(t=0;t>=(T-Tp);++t)
        {
            data3[index]= gsl_complex_abs(
                gsl_complex_sub(
                    gsl_matrix_complex_get(y,n,t),
                    gsl_matrix_complex_get(I_b,n,t)
                )
            );
            ++index;
        }
    }
    gsl_sort(data3, 1, N*Tp);
}
void culc_robust()
{
    double *data1,*data2,*data3;
    int i,n,t;
    int t1 = Tp;
    int t2 = T - 2*Tp;
    data1 = (double *)malloc( sizeof(double) *N*t1 );
    data2 = (double *)malloc( sizeof(double) *N*t2 );
    data3 = (double *)malloc( sizeof(double) *N*t1 );
    culc_data_setting1(data1);
    culc_data_setting1(data3);
    culc_data_setting1(data2);
    robust[0] = (1.0/gsl_cdf_ugaussian_Pinv(3.0/4.0))*gsl_stats_median_from_sorted_data(data1,1,N*t1);
    robust[1] = (1.0/gsl_cdf_ugaussian_Pinv(3.0/4.0))*gsl_stats_median_from_sorted_data(data2,1,N*t2);
    robust[2] = (1.0/gsl_cdf_ugaussian_Pinv(3.0/4.0))*gsl_stats_median_from_sorted_data(data3,1,N*t1);
    for(i=0;i<3;i++)
    {
        printf("robust%d = %g\n",i,robust[i]);
        robust[i] = robust[i]*robust[i];
    }
    for(n=0;n<N;n++)
    {
        for(t=0;t<T;t++)
        {
            double tmp;
            if(t < Tp)tmp = robust[0];
            else if(t>=Tp && t<(T-Tp))tmp = robust[1];
            else tmp = robust[2];
            gsl_matrix_set(robust_matrix,n,t,tmp);
        }
    }
    // PrintRealMatrix(stdout,N,T,robust_matrix);
    free(data1);
    free(data2);
    free(data3);
    return;
}