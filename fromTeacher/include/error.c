#ifndef ERROR_C
#define ERROR_C 1

#include<gsl/gsl_sf_erf.h>

double Qfunction(double x){
	double result;
	result=gsl_sf_erf_Q(x);
	return(result);
}




#endif
