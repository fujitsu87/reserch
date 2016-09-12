#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include "../include/matrix.c"
#include "../include/random_number.c"
#include "../include/blas.c"

//基地局のアンテナの数
#define N 4
//ユーザ数
#define K 4
//離散時間ステップ数
#define T 8
//pilot信号時間長さ
#define Tp 2
//ユーザ分割数
#define DK 2
//pilot信号の値
const double Pk = 1.0;

//ノイズ　10dB
double N0;
//信号が送られる確率
const double rho_k = 1.0;

double a_k;

gsl_matrix_complex* x;
gsl_matrix_complex* h;
gsl_matrix_complex* w;
gsl_matrix_complex* y;

gsl_matrix_complex* z;

gsl_matrix_complex* x_h;
gsl_matrix* x_h_abs2;
gsl_matrix* xi;
gsl_matrix_complex* x_b;
gsl_matrix_complex* xi_b;
gsl_matrix_complex* h_h;
gsl_matrix* h_h_abs2;
gsl_matrix* eta;
gsl_matrix_complex* I_b;
gsl_matrix* zeta;


//--------------------------------初期化関数--------------------------------------

void init_h()
{
	int n,k;
	double sigma =0.5;
	for (n = 0; n < N; ++n)
	{
		for (k = 0; k < K; ++k)
		{
			gsl_complex z = GSLComplexGaussianNoise(sigma);
			// GSL_SET_COMPLEX(&z,1,1);
			gsl_matrix_complex_set(h,n,k,z);
		}
	}
}

void init_x()
{
	int k,t;
	for (k = 0; k < K; ++k)
	{
		for (t = 0; t < T; ++t)
		{

			double tmp[2] = {-1.0,1.0};
			gsl_complex z;
			GSL_SET_COMPLEX(&z,a_k*tmp[UniformBit()],a_k*tmp[UniformBit()]);
			gsl_matrix_complex_set(x ,k,t,z);
		}
	}
}
void init_w()
{
	int n,t;
	for (n = 0; n < N; ++n)
	{
		for (t = 0; t < T; ++t)
		{
			gsl_complex z = GSLComplexGaussianNoise(sqrt(N0));
			// GSL_SET_COMPLEX(&z,1,1);
			gsl_matrix_complex_set(w,n,t,z);
		}
	}
}

void init()
{
	int i,j;
	//ノイズ　10dB
	N0 = Pk/10.0;
	x = gsl_matrix_complex_calloc(K,T);
	h = gsl_matrix_complex_calloc(N,K);
	w = gsl_matrix_complex_calloc(N,T);
	y = gsl_matrix_complex_calloc(N,T);

	z = gsl_matrix_complex_calloc(N,T);

	x_h = gsl_matrix_complex_calloc(K,T);
	x_h_abs2 = gsl_matrix_calloc(K,T);
	xi = gsl_matrix_calloc(K,T);
	x_b = gsl_matrix_complex_calloc(K,T);
	xi_b = gsl_matrix_complex_calloc(K,T);
	h_h = gsl_matrix_complex_calloc(N,K);
	h_h_abs2 = gsl_matrix_calloc(N,K);
	eta = gsl_matrix_calloc(N,K);
	I_b = gsl_matrix_complex_calloc(N,T);
	zeta = gsl_matrix_calloc(N,T);

	a_k = sqrt(Pk/(2*rho_k));
	
	//乱数初期化
	RandomNumberInitialization(0);
	
	//通信路h初期化
	init_h();
	
	//入力値xの初期化
	init_x();
	
	//雑音wの初期化
	init_w();

	//行列計算　y = hx/sqrt(nN) + w
	AB(h,x,y);
	gsl_complex tmp1;
	GSL_SET_COMPLEX(&tmp1,1/sqrt(N),0);
	gsl_matrix_complex_scale(y,tmp1);
	gsl_matrix_complex_add(y,w);

	// init xi
	for (i = 0; i < K; ++i)
	{
		for (j = 0; j < T; ++j)
		{
			gsl_matrix_set(xi ,i,j,1);
		}
	}
	//init eta
	for (i = 0; i < N; ++i)
	{
		for (j = 0; j < K; ++j)
		{
			gsl_matrix_set(eta ,i,j,1);
		}
	}


}

//----------------絶対値の2乗の行列を計算-------------------------------
void abs2_matrix(gsl_matrix_complex *x, gsl_matrix *y,int s1,int s2)
{
	int i,j;
	for (i = 0; i < s1; ++i)
	{
		for (j = 0; j < s2; ++j)
		{
			double tmp = gsl_complex_abs2 (gsl_matrix_complex_get(x,i,j));
			gsl_matrix_set(y ,i,j,tmp);
		}
	}
}
//------------------終了処理-----------------------------------------
void finish()
{
	x = GSLMatrixFree(x); 
	y = GSLMatrixFree(y); 
	h = GSLMatrixFree(h); 
	w = GSLMatrixFree(w); 

	z = GSLMatrixFree(z); 

	x_h = GSLMatrixFree(x_h);
	x_h_abs2 = GSLRealMatrixFree(x_h_abs2); 
	xi = GSLRealMatrixFree(xi);
	x_b = GSLMatrixFree(x_b);
	xi_b = GSLMatrixFree(xi_b);
	h_h = GSLMatrixFree(h_h);
	h_h_abs2 = GSLRealMatrixFree(h_h_abs2);
	eta = GSLRealMatrixFree(eta);
	I_b = GSLMatrixFree(I_b);
	zeta = GSLRealMatrixFree(zeta);
}

//-----------zetaの計算-----------------------------------

void culc_zata()
{
	int n,t,k;
	//culc zata
	gsl_matrix* tmp1;
	gsl_matrix* tmp2;
	gsl_matrix* tmp3;
	tmp1 = gsl_matrix_calloc(N,T);
	tmp2 = gsl_matrix_calloc(N,T);
	tmp3 = gsl_matrix_calloc(N,T);

	abs2_matrix(x_h,x_h_abs2,K,T);
	abs2_matrix(h_h,h_h_abs2,N,K);

	RealAB(eta,xi,tmp1);
	RealAB(eta,x_h_abs2,tmp2);
	RealAB(h_h_abs2,xi,tmp3);
	gsl_matrix_add(zeta,tmp1);
	gsl_matrix_add(zeta,tmp2);
	gsl_matrix_add(zeta,tmp3);
	gsl_matrix_scale(zeta,1/(double)N);

	tmp1 = GSLRealMatrixFree(tmp1);
	tmp2 = GSLRealMatrixFree(tmp2);
	tmp3 = GSLRealMatrixFree(tmp3);

}
//----------------------zの計算----------------------------------------
void culc_z()
{
	gsl_matrix_complex* tmp1;
	tmp1 = gsl_matrix_complex_calloc(N,T);
	gsl_matrix_complex_memcpy(tmp1,y);

	// y - I_b
	gsl_matrix_complex_sub(z,I_b);
	
	//N0 + zeta
	gsl_matrix* tmp2;
	tmp2 = gsl_matrix_calloc(N,T);
	int n,t;
	for (n = 0; n < N; ++n)
	{
		for (t = 0; t < T; ++t)
		{
			double noise = GaussianNoise(sqrt(N0));
			gsl_matrix_set(tmp2,n,t,noise);
		}
	}
	gsl_matrix_add(tmp2,zeta);
	
	//z = (y - I_b)/(N0 + zeta)
	for (n = 0; n < N; ++n)
	{
		for (t = 0; t < T; ++t)
		{
			gsl_complex tmp3 = gsl_complex_div_real(gsl_matrix_complex_get(tmp1,n,t)
											,gsl_matrix_get(tmp2,n,t));
			gsl_matrix_complex_set(z,n,t,tmp3);
		}
	}

	tmp1 = GSLMatrixFree(tmp1);
	tmp2 = GSLRealMatrixFree(tmp2);
}
//------------------通信路推定器----------------------------
void channel_estimation()
{
	culc_zata();
	culc_z();
	
}
//--------------------------------------------------------------
int main()
{
	init();
	// PrintMatrix(stdout,N,T,x_h);
	channel_estimation();
	PrintRealMatrix(stdout,N,T,zeta);
	finish();
	return 0;
}