#include <stdio.h>
#include <math.h>


int main(int argc, char *argv[])
{
    double N,K,T,Tp,NK,TK;
    double tmp1,tmp2;
    printf("input N\n");
    scanf("%lf",&N);
    printf("input K\n");
    scanf("%lf",&K);
    printf("input T\n");
    scanf("%lf",&T);
    NK = N/K;
    TK = T/K;
    tmp1 = NK+TK;
    tmp2 = NK*TK;
    Tp = tmp1/tmp2;
    printf("Tp = %lf\n",Tp*T);
}