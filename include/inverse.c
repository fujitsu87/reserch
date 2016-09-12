#ifndef INVERSE_C
#define INVERSE_C 1

#include<gsl/gsl_matrix.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_permutation.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>
#include"matrix.c"
#include"blas.c"
#include"complex.c"




//#define EPSILON 10e-16
void InverseRealA(gsl_matrix *A,gsl_matrix *B){

	int sigma=0;
	size_t m=A->size1;
	gsl_matrix *LU=gsl_matrix_calloc(m,m);//LU分解したものを格納する行列を確保
	gsl_permutation *p=gsl_permutation_calloc(m);//Aのサイズからpを確保
	
	
	gsl_matrix_memcpy(LU,A);//Aの内容をLUにコピー
	
	gsl_linalg_LU_decomp(LU,p,&sigma);//LU分解
	
	/*ここのコメントアウトをはずして,35行目を消すと複素行列が正則でない場合エラーをだす*/
	/*
	FILE *fp=fopen("process.dat","w");
	double detA=0;
	fprintf(fp,"TransLU:LU\n");
	PrintRealMatrix(fp,m,m,LU);
	detA=gsl_linalg_LU_det(LU,sigma);
	if(detA>EPSILON) gsl_linalg_LU_invert(LU,p,B);
	else if(detA<(-1)*EPSILON) gsl_linalg_LU_invert(LU,p,B);
	else printf("エラー:detAが%+fです\n",detA);
	fclose(fp);
	*/

	gsl_linalg_LU_invert(LU,p,B);
	LU=GSLRealMatrixFree(LU);
	
	gsl_permutation_free(p);
	
}

void InverseA(gsl_matrix_complex *A,gsl_matrix_complex *B){

	int sigma=0;
	size_t m=A->size1;
	gsl_complex detA;
	gsl_matrix_complex *LU=gsl_matrix_complex_calloc(m,m);//LU分解したものを格納する行列を確保
	gsl_permutation *p=gsl_permutation_calloc(m);//Aのサイズからpを確保
	
	
	gsl_matrix_complex_memcpy(LU,A);//Aの内容をLUにコピー
	
	gsl_linalg_complex_LU_decomp(LU,p,&sigma);//LU分解

	/*LU分解の結果を表示したい場合,ここのコメントアウトを外す*/
	/*
	FILE *fp=fopen("process.dat","w");
	fprintf(fp,"TransLU:LU\n");
	PrintRealMatrix(fp,m,m,LU);
	fclose(fp);
	*/
	
	/*行列式が0の場合の例外処理を含む*/
	//detA=gsl_linalg_complex_LU_det(LU,sigma);
	/*
	if(gsl_complex_abs(detA)>EPSILON) gsl_linalg_complex_LU_invert(LU,p,B);
	else if(gsl_complex_abs(detA)<(-1)*EPSILON) gsl_linalg_complex_LU_invert(LU,p,B);
	else printf("エラー:行列Aの固有値が%+f %+fiです\n",detA.dat[0],detA.dat[1]);
	*/

	gsl_linalg_complex_LU_invert(LU,p,B);
	LU=GSLMatrixFree(LU);
	gsl_permutation_free(p);
	
}

void GInverseRealA(gsl_matrix *A,gsl_matrix *B){
	int i=0,j=0;
	size_t n=A->size2;
	
	
	
	gsl_matrix *ATA=gsl_matrix_calloc(n,n);
	gsl_matrix *IATA=gsl_matrix_calloc(n,n);
	
	
	RealATA(A,ATA);
	//gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,A,A,0,ATA);
	InverseRealA(ATA,IATA);
	SymmetricABT(IATA,A,B);

	/*ここのコメントアウトをはずすと,process.datに計算過程を表示する.*/
	/*
	size_t m=A->size1;
	FILE *fp=fopen("process.dat","w");
	gsl_matrix *Recalculation=gsl_matrix_calloc(n,n);

	fprintf(fp,"A\n");
	PrintRealMatrix(fp,m,n,A);
	
	fprintf(fp,"ATA\n");
	PrintRealMatrix(fp,n,n,ATA);

	fprintf(fp,"IATA\n");
	PrintRealMatrix(fp,n,n,IATA);

	RealAB(ATA,IATA,Recalculation);
	fprintf(fp,"ATA*IATA\n");
	PrintRealMatrix(fp,n,n,Recalculation);	
	
	fprintf(fp,"B\n");
	PrintRealMatrix(fp,n,n,B);
	Recalculation=GSLRealMatrixFree(Recalculation);
	fclose(fp);
	*/

	ATA=GSLRealMatrixFree(ATA);
	IATA=GSLRealMatrixFree(IATA);
	
	
	
}

void GInverseA(gsl_matrix_complex *A,gsl_matrix_complex *B){
	int i=0,j=0;

	size_t n=A->size2;
	
	
	
	gsl_matrix_complex *AhA=gsl_matrix_complex_calloc(n,n);
	gsl_matrix_complex *IAhA=gsl_matrix_complex_calloc(n,n);
	
	AHA(A,AhA);
	//gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,UNIT,A,A,ZERO,AhA);
	InverseA(AhA,IAhA);
	HermitianABH(IAhA,A,B);
	
	
	/*ここのコメントアウトをはずすと,process.datに計算過程を表示する.*/
	/*
	size_t m=A->size1;
	FILE *fp=fopen("process.dat","w");
	gsl_matrix_complex *Recalculation=gsl_matrix_complex_calloc(n,n);

	fprintf(fp,"A\n");
	PrintMatrix(fp,m,n,A);
	
	fprintf(fp,"AHA\n");
	PrintMatrix(fp,n,n,AhA);

	fprintf(fp,"IAHA\n");
	PrintMatrix(fp,n,n,IAhA);

	AB(AhA,IAhA,Recalculation);
	fprintf(fp,"AHA*IAHA\n");
	PrintMatrix(fp,n,n,Recalculation);	
	
	fprintf(fp,"B\n");
	PrintMatrix(fp,n,n,B);
	Recalculation=GSLMatrixFree(Recalculation);
	fclose(fp);
	*/

	AhA=GSLMatrixFree(AhA);
	IAhA=GSLMatrixFree(IAhA);
}

#endif
