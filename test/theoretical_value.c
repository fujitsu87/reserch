#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define Pk 1
int main(int argc, char *argv[])
{
    
    if(argc != 4)
    {
        printf("Input error. Please input as below \n" );
        printf("./a.out SNR K Tp\n");
        return 0;
    }

    double N0 = Pk/atof(argv[1]);
    double K =  atof(argv[2]);
    double Tp =  atof(argv[3]);

    // printf("%f %f %f\n",N0,K,Tp);

    double b =1-N0-K/Tp;
    double d = b*b+4*N0;
    double ans[2];
    ans[0] = ((-1)*b+sqrt(d))/2;
    ans[1] = ((-1)*b-sqrt(d))/2;

    double mse = ans[0]/(1+ans[0]);

    printf("sigma^2 = %f %f mse = %f\n",ans[0],ans[1],mse);
    printf("orthogonal_mes = %f\n",N0/(1+N0));
    return 0;
}
