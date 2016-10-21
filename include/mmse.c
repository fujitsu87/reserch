#ifndef MMSE
#define MMSE 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#define EPSILON 1e-13

//関数のパラメータ
struct f_params {
	double snr;	/* SNR */
};

//積分する関数
double f(double x,void *p){
	struct f_params *params = (struct f_params *)p;
	double t1 = 1+tanh(x*params->snr);
	double t2 = x+1;
	double f = t1*t1*exp(-t2*t2*params->snr/2.0);
	return f;
}

double BPSKMMSE(double snr)
{

	//定数1/sqrt(2pisigma^2)
	const double a = 1/sqrt(2*M_PI/snr);

	//計算結果と推定誤差
	double result, error;

	//数値積分の分割数
	const size_t n = 100000; 
	//パラメータの代入
	struct f_params params =
	{
		snr
	};
	//gsl　functionの定義
	gsl_function F;
	F.function = &f;
	F.params = &params;

	// n個の区間での積分結果と推定誤差を倍精度で保持するための作業領域を確保する。
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(n);

	//無限区間積分
	gsl_integration_qagi(&F,0, EPSILON, n, w, &result, &error);
	
	//領域の解放
	gsl_integration_workspace_free(w);
	return a*result;
}
#endif
