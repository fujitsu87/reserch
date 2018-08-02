#ifndef _MEDIAN_C 
#define _MEDIAN_C 1 

#include<gsl/gsl_sort.h>
#include<gsl/gsl_statistics_double.h>

#include"../../../include/vector.c"
#include"../../../include/random_number.c"

double Median(int n,gsl_vector *x)
{
  int i;
  double *x_tmp,med; 

  x_tmp = DVectorAlloc(0,n-1);
  for(i=0;i<n;i++) {
    x_tmp[i] = gsl_vector_get(x,i);
  }

  gsl_sort(x_tmp,1,n);
  med = gsl_stats_median_from_sorted_data(x_tmp,1,n);
  x_tmp = DVectorFree(x_tmp,0,n-1);
  return med; 
}

void Swap(double *a,double *b)
{
  double tmp;

  tmp = a[0];
  //printf("swap %g %g\n",a[0],b[0]);
  a[0] = b[0];
  //printf("swap %g %g\n",a[0],b[0]);
  b[0] = tmp;
  //printf("swap %g %g\n",a[0],b[0]);
}

int Pivot(int n_first,int n_end)
{
  int pivot;

  pivot = UniformInt(n_end-n_first+1) + n_first;
  return pivot;
}

int Partition(double *x,int n_first,int n_end,int pivot)
{
  int n,index;
  double x_pivot;

  x_pivot = x[pivot];
  Swap(&x[pivot],&x[n_end]);

  index = n_first;
  for(n=n_first;n<n_end;n++) {
    if(x[n]<=x_pivot) {
      Swap(&x[index],&x[n]);
      index++;
    }
  }
  Swap(&x[n_end],&x[index]);

  return index;
}

void Selection_End(double *x,int pivot,int n_first)
{
  int index,n;
  double x_max; 

  x_max = x[n_first]; index = n_first;
  for(n=n_first+1;n<pivot;n++) {
    if(x_max<x[n]) {
      x_max = x[n]; index = n;
    }
  }
  Swap(&x[pivot-1],&x[index]);
}

double Selection(double *x,int k,int first,int end)
{
  int pivot,pivot_new,n,n_first,n_first_pre,n_end,n_end_pre;

  n_first = n_first_pre = first; 
  n_end = end;
  for(;;) {
    pivot = Pivot(n_first,n_end);
    //printf("pivot=%g first=%d end=%d\n",x[pivot],n_first,n_end);
    pivot_new = Partition(x,n_first,n_end,pivot);

    /*
    for(n=first;n<=end;n++) {
      printf("%g ",x[n]);
    }
    printf("\n");
    */

    //printf("n_first_pre=%d\n",n_first_pre);
    if(k==pivot_new) {
      if((end-first+1)%2==0) {

	//printf("%d %d\n",n_first_pre,pivot_new);
	Selection_End(x,pivot_new,n_first_pre);

	/*
	printf("Termination\n");
	for(n=first;n<=end;n++) {
	  printf("%g ",x[n]);
	}
	printf("\n");
	*/

	return (x[pivot_new-1]+x[pivot_new])*0.5;
      }
      else return x[pivot_new];
    }
    else {
      if(k<pivot_new) {
	n_end = pivot_new - 1;
      }
      else {
	n_first = pivot_new + 1;
	if(k!=n_first) n_first_pre = n_first; 
      }
    }
  }
}

double QMedian(int n,gsl_vector *x_vec)
{
  int i;
  double *x,median;

  x = DVectorAlloc(1,n);

  for(i=0;i<n;i++) {
    x[i] = gsl_vector_get(x_vec,i);
    //printf("%g ",x[i]);
  }
  //printf("\n");

  median = Selection(x,n/2,0,n-1);

  /*
  for(i=0;i<n;i++) {
    printf("%g ",x[i]);
  }
  printf("\n");
  printf("%g\n",median);
  */

  //printf("%g\n",median-Median(n,x_vec));

  x = DVectorFree(x,1,n);
}


#endif
