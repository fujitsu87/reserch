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
			if(gsl_matrix_get(A,i,j)>=0) fprintf(fp,"+%g   \t",gsl_matrix_get(A,i,j));
			else fprintf(fp,"-%g   \t",-1*gsl_matrix_get(A,i,j));//fpで示されたファイルに行列Aを出力
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
	 		if(GSL_REAL(z)>=0&&GSL_IMAG(z)>=0) fprintf(fp,"+%g\t+i%g\t\t",GSL_REAL(z),GSL_IMAG(z));//実数成分が正,虚数成分が正の場合
			else if(GSL_REAL(z)<=0&&GSL_IMAG(z)>=0)fprintf(fp,"-%g\t+i%g\t\t",-1*GSL_REAL(z),GSL_IMAG(z));//実数成分が負,虚数成分が正の場合
			else if(GSL_REAL(z)>=0&&GSL_IMAG(z)<=0)fprintf(fp,"+%g\t-i%g\t\t",GSL_REAL(z),-1*GSL_IMAG(z));//実数成分が正,虚数成分が負の場合
			else fprintf(fp,"-%g\t-i%g\t\t",-1*GSL_REAL(z),-1*GSL_IMAG(z));//両方が負
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");//1行空白をいれる

}

int **IMatrix(unsigned int m,unsigned int n,int **mat){
  int i,j;
  mat = (int**)calloc(m,sizeof(int *) );
  mat[0]=(int*)calloc(m*n,sizeof(int));
  for (i=1;i<m;i++) {
     mat[i] = mat[0]+i*n;	//整数型でのmxnの配列を動的確保
  }

  return mat;
}
double **IMatrixf(unsigned int m,unsigned int n,double **mat){
  int i,j;
  mat = (double**)calloc(m,sizeof(double *) );
  mat[0]=(double*)calloc(m*n,sizeof(double));
  for (i=1;i<m;i++) {
     mat[i] = mat[0]+i*n;	//浮動小数点型でのmxnの配列を動的確保
  }

  return mat;
}
void IMatrixFree(int **mat){	//整数型で確保された配列のメモリ領域を解放
  free(mat[0]);
  free(mat);
}
void IMatrixfFree(double **mat){//浮動小数点型で確保された配列のメモリ領域を解放
  free(mat[0]);
  free(mat);
}

#endif
