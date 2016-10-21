#ifndef SORT_C
#define SORT_C 1

#include <stdlib.h>

/* ソート関数 */
void DSortDes(unsigned int n, double array[n]);

/* 比較関数(static) */
static int DCompareDes(const void *a, const void *b);


/* ソート関数本体 */
//----------------double型配列の降順ソート----------------
void DSortDes(unsigned int n, double array[n])
{
	qsort(array, n, sizeof(double), DCompareDes);
}

/* 比較関数(static)本体 */
//----------------double型配列の降順ソート用比較関数----------------
static int DCompareDes(const void *a, const void *b)
{
	if(*(double*)b < *(double*)a)
		return -1;
	else if(*(double*)b > *(double*)a)
		return 1;
	else
		return 0;
}

#endif