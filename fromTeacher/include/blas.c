#ifndef BLAS_C
#define BLAS_C 1

#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>
#include<gsl/gsl_matrix.h>
#include"complex.c"
#include"matrix.c"

/*実ベクトルのユークリッドノルムを返す関数*/
double RealNorm2(gsl_vector *x){
	return(gsl_blas_dnrm2(x));
}

/*複素ベクトルのユークリッドノルムを返す関数*/
double Norm2(gsl_vector_complex *x){
	return(gsl_blas_dznrm2(x));
}

/*2つの実ベクトルの標準内積を返す関数*/
double RealxTy(gsl_vector *x,gsl_vector *y){
	double result;//ポインタで宣言しなければいけないことに注意
	gsl_blas_ddot(x,y,&result);
	return(result);
}

/*2つの複素ベクトルの標準内積を返す関数*/
gsl_complex xHy(gsl_vector_complex *x,gsl_vector_complex *y){
	gsl_complex result;
	gsl_blas_zdotc(x,y,&result);
	return(result);
}

/*実数空間においてy=Axを計算する関数*/
void RealAx(gsl_matrix *A,gsl_vector *x,gsl_vector *y){
	gsl_blas_dgemv(CblasNoTrans,1.0,A,x,0,y);
}

/*実数空間においてy=A^(T)xを計算する関数*/
void RealATx(gsl_matrix *A,gsl_vector *x,gsl_vector *y){
	gsl_blas_dgemv(CblasTrans,1.0,A,x,0,y);
}

/*複素空間においてy=Axを計算する関数*/
void Ax(gsl_matrix_complex *A,gsl_vector_complex *x,gsl_vector_complex *y){
	gsl_complex alpha;
	gsl_complex beta;
	alpha=UNIT;
	beta=ZERO;

	gsl_blas_zgemv(CblasNoTrans,alpha,A,x,beta,y);
}

/*複素空間においてy=A^(H)xを計算する関数*/
void AHx(gsl_matrix_complex *A,gsl_vector_complex *x,gsl_vector_complex *y){
	gsl_complex alpha;
	gsl_complex beta;
	alpha=UNIT;
	beta=ZERO;

	gsl_blas_zgemv(CblasConjTrans,alpha,A,x,beta,y);
}

/*実空間において対象行列Aのy=Axを計算する関数*/
void SymmetricAx(gsl_matrix *A,gsl_vector *x,gsl_vector *y){
	gsl_blas_dsymv(CblasLower,1.0,A,x,0,y);
}

/*エルミート行列Aでy=Axを計算する関数*/
void HermitianAx(gsl_matrix_complex *A,gsl_vector_complex *x,gsl_vector_complex *y){
	gsl_complex alpha;
	gsl_complex beta;
	alpha=UNIT;
	beta=ZERO;

	gsl_blas_zhemv(CblasLower,alpha,A,x,beta,y);
}

/*実空間において実行列の積C=ABを計算する関数*/
void RealAB(gsl_matrix *A,gsl_matrix *B,gsl_matrix *C){
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,A,B,0,C);
}

/*実空間において実行列の積C=A^(T)Bを計算する関数*/
void RealATB(gsl_matrix *A,gsl_matrix *B,gsl_matrix *C){
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,A,B,0,C);
}

/*実空間において実行列の積C=AB^(T)を計算する関数*/
void RealABT(gsl_matrix *A,gsl_matrix *B,gsl_matrix *C){
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,A,B,0,C);
}

/*複素空間において複素行列の積C=ABを計算する関数*/
void AB(gsl_matrix_complex *A,gsl_matrix_complex *B,gsl_matrix_complex *C){
	gsl_complex alpha;
	gsl_complex beta;
	alpha=UNIT;
	beta=ZERO;
	gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,alpha,A,B,beta,C);
}

/*複素空間において複素行列の積C=A^(T)Bを計算する関数*/
void AHB(gsl_matrix_complex *A,gsl_matrix_complex *B,gsl_matrix_complex *C){
	gsl_complex alpha;
	gsl_complex beta;
	alpha=UNIT;
	beta=ZERO;
	gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,alpha,A,B,beta,C);
}

/*複素空間において複素行列の積C=AB^(T)を計算する関数*/
void ABH(gsl_matrix_complex *A,gsl_matrix_complex *B,gsl_matrix_complex *C){
	gsl_complex alpha;
	gsl_complex beta;
	alpha=UNIT;
	beta=ZERO;
	gsl_blas_zgemm(CblasNoTrans,CblasConjTrans,alpha,A,B,beta,C);
}

/*実空間においてB=AA^(T)を計算する関数*/
void RealAAT(gsl_matrix *A,gsl_matrix *B){
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,A,A,0,B);
}

/*実空間においてB=A^(T)Aを計算する関数*/
void RealATA(gsl_matrix *A,gsl_matrix *B){
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,A,A,0,B);
}

/*複素空間においてB=AA^(H)を計算する関数*/
void AAH(gsl_matrix_complex *A,gsl_matrix_complex *B){
	gsl_complex alpha;
	gsl_complex beta;
	alpha=UNIT;
	beta=ZERO;
	gsl_blas_zgemm(CblasNoTrans,CblasConjTrans,alpha,A,A,beta,B);
}

/*複素空間においてB=A^(H)Aを計算する関数*/
void AHA(gsl_matrix_complex *A,gsl_matrix_complex *B){
	gsl_complex alpha;
	gsl_complex beta;
	alpha=UNIT;
	beta=ZERO;
	gsl_blas_zgemm(CblasConjTrans,CblasNoTrans,alpha,A,A,beta,B);
}

/*実空間において対象行列Aのy=AB^(T)を計算する関数*/
void SymmetricABT(gsl_matrix *A,gsl_matrix *B,gsl_matrix *C){
	gsl_matrix *BT=gsl_matrix_calloc(B->size2,B->size1);
	gsl_matrix_transpose_memcpy(BT,B);
	gsl_blas_dsymm(CblasLeft,CblasLower,1.0,A,BT,0,C);
	BT=GSLRealMatrixFree(BT);
}

/*複素空間においてエルミート行列AのC=AB^(H)を計算する関数*/
void HermitianABH(gsl_matrix_complex *A,gsl_matrix_complex *B,gsl_matrix_complex *C){
	int i,j;
	gsl_complex z;
	gsl_matrix_complex *BH=gsl_matrix_complex_calloc(B->size2,B->size1);
	//FILE *fp=fopen("process.dat","w");
	ConjugateTranspose(B,BH);
	//PrintMatrix(fp,BH->size1,BH->size2,BH);
	//AB(A,BH,C);
	gsl_blas_zhemm(CblasLeft,CblasLower,UNIT,A,BH,ZERO,C);
	BH=GSLMatrixFree(BH);
	//fclose(fp);
}


#endif
