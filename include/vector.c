#ifndef VECTOR_C
#define VECTOR_C 1

#include <gsl/gsl_vector.h>
#include <stdio.h>
#include <gsl/gsl_complex.h>

//----------------実ベクトルの解放----------------
gsl_vector *GSLRealVectorFree(gsl_vector *v)
{
	if(v != NULL)
		gsl_vector_free(v);
	return NULL;
}

//----------------複素ベクトルの解放----------------
gsl_vector_complex *GSLVectorFree(gsl_vector_complex *v)
{
	if(v != NULL)
		gsl_vector_complex_free(v);
	return NULL;
}

//----------------実ベクトルのファイル出力----------------
void PrintRealVector(FILE* fp, unsigned int n, gsl_vector *v)
{
	int i;
	int N = (n < v->size) ? n : v->size;

	for(i = 0; i < N; i++)
	{
		fprintf(fp,"%+f",gsl_vector_get(v,i));
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
}

//----------------複素ベクトルのファイル出力----------------
void PrintVector(FILE* fp, unsigned int n, gsl_vector_complex *v)
{
	int i;
	int N = (n < v->size) ? n : v->size;
	gsl_complex z;

	for(i = 0; i < N; i++)
	{
		z = gsl_vector_complex_get(v,i);
		fprintf(fp,"%+f",GSL_REAL(z));		//実部
		fprintf(fp,"\t%+f i",GSL_IMAG(z));	//虚部
		fprintf(fp,"\n");					//改行
	}
	fprintf(fp,"\n");
}

#endif