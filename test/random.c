#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "../include/matrix.c"
#include "../include/random_number.c"
#include "../include/blas.c"
#include "../include/inverse.c"


int main(int argc, char *argv[])
{
    RandomNumberInitialization(2);
    printf("%g\n",gsl_rng_uniform(RAN));
    printf("%g\n",gsl_rng_uniform(RAN));


    return 0;
}