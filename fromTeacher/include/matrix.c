#ifndef MATRIX_C
#define MATRIX_C 1
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_complex.h>

/*実行列の開放を行う関数*/
gsl_matrix *GSLRealMatrixFree(gsl_matrix*A)
{
	if(A==NULL){
		//printf("NULLを返しましたよ\n");
		return(NULL);
		}
	else{
		//printf("行列を開放しましたよ\n");
		gsl_matrix_free(A);//与えられた行列を開放する
		return(NULL);
		}
}

/*実行列をテキストファイルに出力する関数*/
gsl_matrix *PrintRealMatrix(FILE *fp,unsigned int m,unsigned int n,gsl_matrix *A)
{	
	unsigned int i,j;
	

	for(i=0;i<m;i++){
		for(j=0;j<n;j++){
			if(gsl_matrix_get(A,i,j)>=0) fprintf(fp,"%+f   \t",gsl_matrix_get(A,i,j));
			else fprintf(fp,"%+f   \t",gsl_matrix_get(A,i,j));//fpで示されたファイルに行列Aを出力
		}
		fprintf(fp,"\n");
	}		
	fprintf(fp,"\n");//最後の1行をあける
}

/*虚数行列の開放を行う関数*/
gsl_matrix_complex *GSLMatrixFree(gsl_matrix_complex *A)
{
		if(A==NULL){
		//printf("NULLを返しましたよ\n");
		return(NULL);
		}
	else{
		gsl_matrix_complex_free(A);//与えられた行列を開放する
		return(NULL);
		//printf("行列を開放しましたよ\n");
		}
	
}

/*虚数行列をテキストファイルに出力する関数*/
gsl_matrix *PrintMatrix(FILE *fp,unsigned int m,unsigned int n,gsl_matrix_complex *A)
{	
	unsigned int i,j;
	gsl_complex z;//虚数を定義

	for(i=0;i<m;i++){
		for(j=0;j<n;j++){
			z=gsl_matrix_complex_get(A,i,j);
			/*
	 		if(GSL_REAL(z)>=0&&GSL_IMAG(z)>=0) fprintf(fp,"+%g\t+i%g\t\t",GSL_REAL(z),GSL_IMAG(z));//実数成分が正,虚数成分が正の場合
			else if(GSL_REAL(z)<=0&&GSL_IMAG(z)>=0)fprintf(fp,"-%g\t+i%g\t\t",-1*GSL_REAL(z),GSL_IMAG(z));//実数成分が負,虚数成分が正の場合
			else if(GSL_REAL(z)>=0&&GSL_IMAG(z)<=0)fprintf(fp,"+%g\t-i%g\t\t",GSL_REAL(z),-1*GSL_IMAG(z));//実数成分が正,虚数成分が負の場合
			else fprintf(fp,"-%g\t-i%g\t\t",-1*GSL_REAL(z),-1*GSL_IMAG(z));//両方が負
			*/
			fprintf(fp,"%+f\t%+fi\t",GSL_REAL(z),GSL_IMAG(z));
		}
		fprintf(fp,"\n");
	}	
	fprintf(fp,"\n");//1行空白をいれる	
	
}
gsl_complex gsl_complex_conjugate (gsl_complex z);
void ConjugateTranspose(gsl_matrix_complex *A,gsl_matrix_complex *AH)
{
	int i,j;
	gsl_complex z;
	gsl_matrix_complex_transpose_memcpy(AH,A);
	for(i=0;i<AH->size1;i++){
		for(j=0;j<AH->size2;j++){
			z = gsl_matrix_complex_get(AH,i,j);
			z = gsl_complex_conjugate(z);
			gsl_matrix_complex_set(AH,i,j,z);
		}
	}
}

#endif