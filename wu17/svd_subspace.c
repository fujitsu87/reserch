#include "../include/svd.c"
#include "../include/matrix.c"
#include "../include/complex.c"

void get_s(gsl_matrix_complex *U)
{
    int k,n;
    for(n = 0; n < N; ++n)
        for(k = 0; k < K/DK; ++k)
            gsl_matrix_complex_set (s, n, k,gsl_matrix_complex_get (U, n, k));
    return;
} 
void svd()
{
    gsl_matrix_complex *U = gsl_matrix_complex_alloc(N, N);

    double sigma[N];
    unsigned int rank;
    int i;

    //Using functions
    SVDLeft(y, sigma, U);
    get_s(U);
    //Outputing
    // PrintMatrix(stdout, N, T, y);
    // fprintf(stdout, "---␣SingularValues␣of␣A␣(using␣’SingularValue’)␣---\n");
    // for(i = 0; i < T; i++)
    // fprintf(stdout, "%+.18f\t", sigma[i]);
    // fprintf(stdout, "\n\n---␣U␣,␣LeftSingularVectors␣of␣A␣(using␣’SVDLeft’)␣---\n");
    // PrintMatrix(stdout, N, N, U);
    // fprintf(stdout, "\n\n---S---\n");
    // PrintMatrix(stdout, N, K, s);
    //Freeing

    // gsl_matrix_complex *Uy = gsl_matrix_complex_alloc(N, T);
    // AHB(U,y,Uy);
    // fprintf(stdout, "\n\n---Uy---\n");
    // PrintMatrix(stdout, N, T, Uy);
    // Uy = GSLMatrixFree(Uy);
    
    // gsl_matrix_complex *xxh = gsl_matrix_complex_alloc(K, K);
    // ABH(x_h,x_h,xxh);
    // fprintf(stdout, "\n\n---xxh---\n");
    // PrintMatrix(stdout, K, K, xxh);
    // xxh = GSLMatrixFree(xxh);

    U = GSLMatrixFree(U);
}