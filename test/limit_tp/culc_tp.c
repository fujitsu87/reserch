#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    if(argc != 4)
    {
      printf("Input error. please input as below\n");
      printf("./a.out N K T\n");
      return 0;
    }
    double N,K,T,Tp,NK,TK;
    double tmp1,tmp2;
    //printf("input N\n");
    //scanf("%lf",&N);
    N = atof(argv[1]);
    //printf("input K\n");
    //scanf("%lf",&K);
    K = atof(argv[2]);
    //printf("input T\n");
    //scanf("%lf",&T); 
    T = atof(argv[3]);
    int i;
    for(i=16;i<=2048;i=i*2){
	T = (double)i;
    	NK = N/K;
    	TK = T/K;
    	tmp1 = NK+TK;
    	tmp2 = NK*TK;
    	Tp = tmp1/tmp2;
    	printf("%lf %lf %lf\n",T,Tp*T,Tp);
    }
}
