#ifndef RANDOM_NUMBER
#define RANDOM_NUMBER 1

#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_complex.h>
#include <complex.h>
#include <stdlib.h>


//グローバル変数
gsl_rng *RAN;

//----------------乱数生成器を初期化-------------------------
void RandomNumberInitialization(unsigned int SEED )
{
	gsl_rng_env_setup();

	RAN=gsl_rng_alloc(gsl_rng_default);

	if(SEED==0){
		gsl_rng_set(RAN,time(NULL));   //歴時間を種として設定
	}
	else {
		gsl_rng_set(RAN,SEED);        //与えられた自然数を種に設定
	}
}

//--------------------二元一様乱数を生成----------------------
unsigned char UniformBit()
{
	unsigned char unibin;
	unibin=gsl_rng_uniform_int(RAN,2);
	return unibin;
}

//-------------------実ガウスノイズを生成----------------------
double GaussianNoise(double sigma)
{
	double gaussn;
	gaussn=gsl_ran_gaussian(RAN,sigma);   //標準偏差sigmaのガウス雑音を生成
	return gaussn;
}

//------------------複素ガウスノイズを生成------------------------

complex ComplexGaussianNoise(double sigma)
{
	complex z;//z=x+iyを定義
	double standard_dev = sigma*M_SQRT1_2;
	z=gsl_ran_gaussian(RAN,standard_dev)+I*gsl_ran_gaussian(RAN,standard_dev);//平均0,sigma/root2のガウス雑音を代入
	return z;
}

//-------------------GSL用の複素ガウスノイズを生成----------------
gsl_complex GSLComplexGaussianNoise(double sigma)
{
	gsl_complex z;
	double standard_dev = sigma*M_SQRT1_2;
	double y = gsl_ran_gaussian(RAN,standard_dev);
	double x = gsl_ran_gaussian(RAN,standard_dev);
	GSL_SET_COMPLEX(&z, x, y);
	return z;
}


#endif