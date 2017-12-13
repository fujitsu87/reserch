#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    if(argc != 3)
    {
      printf("Input error. please input as below\n");
      printf("./a.out N T\n");
      return 0;
    }
    double N,T,NT;
    double ans;
    N = atof(argv[1]);
    T = atof(argv[2]);

    NT =N*T;
    ans = NT/(N+T);

    printf("%lf\n",ans);
}
