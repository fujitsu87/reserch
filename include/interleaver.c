#include "../include/random_number.c"

extern gsl_rng *RAN;

void Permutation(int n,int P[n]){
  int i;

  for(i=0;i<n;i++){
    P[i]=i;
  }
  gsl_ran_shuffle(RAN,P,n,sizeof(int));

}

void RealInterleaver(int n,double x[n],double y[n],int P[n]){
  int i;
  for(i=0;i<n;i++){
    y[P[i]]=x[i];
  }
}
