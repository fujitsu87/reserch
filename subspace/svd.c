#include "../include/svd.c"
#include "../include/matrix.c"
#include "../include/complex.c"

const size_t M = 3; //number of rows
const size_t N = 2; //number of columns
const size_t K = 2; //min(M, N)

int main(void)
{
    gsl_matrix_complex * A = gsl_matrix_complex_alloc(M, N);
    gsl_matrix_complex * Ul = gsl_matrix_complex_alloc(M, M);
    gsl_matrix_complex * Vr = gsl_matrix_complex_alloc(N, N);
    gsl_matrix_complex * U = gsl_matrix_complex_alloc(M, M);
    gsl_matrix_complex * V = gsl_matrix_complex_alloc(N, N);
    gsl_matrix_complex * USVH = gsl_matrix_complex_alloc(M, N);
    double sigma[K] ,sigma_L[M], sigma_R[N], sigma_SVD[K];
    unsigned int rank;
    //Making A
    int i,j;
    for(i = 0; i < M; i++){
        for(j = 0; j < N; j++){
        gsl_matrix_complex_set(A, i, j, gsl_complex_rect(1-i+j, 2+i-j));
        }
    }

    //Using functions
    SingularValue(A, sigma);
    rank = SVDRank(A);
    SVDLeft(A, sigma_L, Ul);
    SVDRight(A, sigma_R, Vr);
    SVD(A, sigma_SVD, U, V);

    //Making USVH
    gsl_matrix_complex * S = gsl_matrix_complex_calloc(M, N);
    gsl_matrix_complex * US = gsl_matrix_complex_alloc(M, N);
    for(i = 0; i < K; i++)
    gsl_matrix_complex_set(S, i, i, gsl_complex_rect(sigma_SVD[i], 0));
    AB(U, S, US);
    ABH(US, V, USVH);
    S = GSLMatrixFree(S);
    US = GSLMatrixFree(US);
    //Outputing
    fprintf(stdout, "---␣Input␣Matrix␣A␣(rank␣=␣%d)␣---\n", rank);

    PrintMatrix(stdout, M, N, A);
    fprintf(stdout, "---␣SingularValues␣of␣A␣(using␣’SingularValue’)␣---\n");
    for(i = 0; i < K; i++)
    fprintf(stdout, "%+.18f\t", sigma[i]);
    fprintf(stdout, "\n\n---␣SingularValues␣of␣A␣(using␣’SVDLeft’)␣---\n");
    for(i = 0; i < M; i++)
    fprintf(stdout, "%+.18lf\t", sigma_L[i]);
    fprintf(stdout, "\n\n---␣SingularValues␣of␣A␣(using␣’SVDRight’)␣---\n");
    for(i = 0; i < N; i++)
    fprintf(stdout, "%+.18lf\t", sigma_R[i]);
    fprintf(stdout, "\n\n---␣SingularValues␣of␣A␣(using␣’SVD’)␣---\n");
    for(i = 0; i < K; i++)
    fprintf(stdout, "%+.18lf\t", sigma_SVD[i]);

    fprintf(stdout, "\n\n---␣U␣,␣LeftSingularVectors␣of␣A␣(using␣’SVDLeft’)␣---\n");
    PrintMatrix(stdout, M, M, Ul);
    fprintf(stdout, "---␣V␣,␣RightSingularVectors␣of␣A␣(using␣’SVDRight’)␣---\n");

    PrintMatrix(stdout, N, N, Vr);
    fprintf(stdout, "---␣U␣,␣LeftSingularVectors␣of␣A␣(using␣’SVD’)␣---\n");
    PrintMatrix(stdout, M, M, U);
    fprintf(stdout, "---␣V␣,␣RightSingularVectors␣of␣A␣(using␣’SVD’)␣---\n");
    PrintMatrix(stdout, N, N, V);

    fprintf(stdout, "---␣USVH␣---\n");
    PrintMatrix(stdout, M, N, USVH);
    //Freeing
    A = GSLMatrixFree(A);
    Ul = GSLMatrixFree(Ul);
    Vr = GSLMatrixFree(Vr);
    U = GSLMatrixFree(U);
    V = GSLMatrixFree(V);
    USVH = GSLMatrixFree(USVH);
}