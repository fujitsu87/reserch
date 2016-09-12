#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include "../include/matrix.c"
#include "../include/random_number.c"
#include "../include/blas.c"

//基地局のアンテナの数
#define N 30
//ユーザ数
#define K 50
//離散時間ステップ数
#define T 150
//pilot信号時間長さ
#define Tp 50
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
gsl_matrix* xi_b;
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

			//左下のpilot信号
			if(k >= K/DK && t < Tp)
			{
				gsl_matrix_complex_set(x_h,k,t,z);
			}
						//左下のpilot信号
			else if(k < K/DK && t >= T-Tp)
			{
				gsl_matrix_complex_set(x_h,k,t,z);
			}
		
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
	xi_b = gsl_matrix_calloc(K,T);
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
	
	//入力値xの初期化 pilot信号をx_hにも代入
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
			gsl_matrix_set(xi_b ,i,j,1);
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
	xi_b = GSLRealMatrixFree(xi_b);
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
//------------------A関数-------------------------------------------
double diff_fk(double u,double v)
{
	double ans;
	double e1,e2,e3,t1,t2;
	e1 = exp(2*a_k*u/v);
	e2 = exp(-2*a_k*u/v);
	e3 = exp(a_k*a_k/v);

	t1 = 2*a_k/v*(e1*(a_k-1)+e2*(a_k-1));
	t2 = e1 + e2 + 2 * (1/rho_k - 1) * e3;

	t2 = t2 * t2;

	ans = t1/t2; 

	return ans;
}

gsl_complex A(gsl_complex arg,int k,int t)
{
	gsl_complex ans;
	double v = gsl_matrix_get(xi_b,k,t);
	ans.dat[0] = GSL_REAL(arg) 
				* diff_fk(GSL_REAL(gsl_matrix_complex_get(x_h,k,t))
							,v
					);
	ans.dat[1] = GSL_IMAG(arg) 
				* diff_fk(GSL_IMAG(gsl_matrix_complex_get(x_h,k,t))
							,v
					);
	return ans;
}

//----------------I_bの計算--------------------------------
void culc_I()
{
	int n,t,k;
	gsl_complex tmp;
	gsl_complex tmp2;
	gsl_complex tmp3;
	gsl_complex tmp4;
	abs2_matrix(x_h,x_h_abs2,K,T);

	gsl_matrix_complex *h_h_conju;
	h_h_conju = gsl_matrix_complex_calloc(K,N);
	ConjugateTranspose(h_h,h_h_conju);

	// PrintMatrix(stdout,K,N,h_h_conju);

	// 1/sqrt(N) h_h * x_h
	AB(h_h,x_h,I_b);
	double n_sqrt = sqrt((double)N);
	GSL_SET_COMPLEX(&tmp,1/n_sqrt,0);
	gsl_matrix_complex_scale(I_b,tmp);

	//xi_b * A()
	for (n = 0; n < N; ++n)
	{
		for (t = 0; t < T; ++t)
		{
			GSL_SET_COMPLEX(&tmp2,0,0);
			for (k = 0; k < K; ++k)
			{
				GSL_SET_COMPLEX(&tmp3,0,0);
				//xiを取り出す
				GSL_SET_COMPLEX(&tmp3,gsl_matrix_get(xi_b,k,t),0);
				
				//Aktの計算
				tmp4 = A(
							gsl_complex_mul(
								gsl_matrix_complex_get(h_h_conju,k,n),
								gsl_matrix_complex_get(z,n,t)
							)
							,k
							,t
						);
				//xi_b * A
				tmp3 = gsl_complex_mul(tmp3,tmp4);
				//h_hとかけ合わせる
				tmp3 = gsl_complex_mul(gsl_matrix_complex_get(h_h,n,k),tmp3);

				//eta * x_h^2 * znt
				tmp4 = gsl_complex_mul_real(
							gsl_matrix_complex_get(z,n,t)
							,gsl_matrix_get(eta,n,k)*gsl_matrix_get(x_h_abs2,k,t)
						);


				tmp2.dat[0] += tmp3.dat[0] + tmp4.dat[0];
				tmp2.dat[1] += tmp3.dat[1] + tmp4.dat[1];
			}

			tmp2.dat[0] = tmp2.dat[0]/(double)N;
			tmp2.dat[1] = tmp2.dat[1]/(double)N;

			//I_bから引く
			tmp3 = gsl_complex_sub(gsl_matrix_complex_get(I_b,n,t),tmp2);
			//I_bにset
			gsl_matrix_complex_set(I_b,n,t,tmp3);

		}
	}

	h_h_conju = GSLMatrixFree(h_h_conju);
}
//------------------culc_h_h-----------------------------------
void culc_h_h()
{
	int n,t,k;
	gsl_complex tmp1;
	gsl_complex tmp2;
	gsl_complex tmp3;
	gsl_complex sec1;
	gsl_complex sec2;
	gsl_complex sec3;
	GSL_SET_COMPLEX(&sec1,0,0);
	GSL_SET_COMPLEX(&sec2,0,0);
	GSL_SET_COMPLEX(&sec3,0,0);


	gsl_matrix_complex *x_h_conju;
	x_h_conju = gsl_matrix_complex_calloc(T,K);
	ConjugateTranspose(x_h,x_h_conju);

	gsl_matrix_complex *h_h_conju;
	h_h_conju = gsl_matrix_complex_calloc(K,N);
	ConjugateTranspose(h_h,h_h_conju);


	double n_sqrt = sqrt((double)N);


	for (n = 0; n < N; ++n)
	{
		for (k = 0; k < K; ++k)
		{
			GSL_SET_COMPLEX(&sec1,0,0);
			GSL_SET_COMPLEX(&sec2,0,0);
			GSL_SET_COMPLEX(&sec3,0,0);
			for (t = 0; t < T; ++t)
			{
				GSL_SET_COMPLEX(&tmp1,0,0);
				
				//h_h * x_h_conju　第1項の計算
				
				tmp1 = gsl_complex_mul(
							gsl_matrix_complex_get(x_h_conju,k,n)
							,gsl_matrix_complex_get(z,n,t)
						);

				sec1.dat[0] += tmp1.dat[0];
				sec1.dat[1] += tmp1.dat[1];
				// printf("%f %f\n",sec1.dat[0],sec1.dat[1]);
				//第3項の計算
				GSL_SET_COMPLEX(&tmp2,0,0);
				//Aktの計算
				tmp2 = A(
							gsl_complex_mul(
								gsl_matrix_complex_get(h_h_conju,k,n),
								gsl_matrix_complex_get(z,n,t)
							)
							,k
							,t
						);
				tmp2 = gsl_complex_conjugate(tmp2);
				tmp2 = gsl_complex_mul_real(
							tmp2
							,gsl_matrix_get(xi_b,k,t)
						);
				tmp2 = gsl_complex_mul(
							gsl_matrix_complex_get(z,n,t)
							,tmp2
						);
				sec3.dat[0] += tmp2.dat[0];
				sec3.dat[1] += tmp2.dat[1];

			}
			sec1.dat[0] = sec1.dat[0]/n_sqrt;
			sec1.dat[1] = sec1.dat[1]/n_sqrt;
			sec1 = gsl_complex_mul_real(sec1,gsl_matrix_get(eta,n,k));

			sec2 = gsl_complex_mul_real(
					gsl_matrix_complex_get(h_h,n,k)
					,(1-gsl_matrix_get(eta,n,k))
				);

			sec3.dat[0] = sec3.dat[0]/(double)N;
			sec3.dat[1] = sec3.dat[1]/(double)N;
			sec3 = gsl_complex_mul_real(sec3,-gsl_matrix_get(eta,n,k));

			GSL_SET_COMPLEX(&tmp3,0,0);
			tmp3 = gsl_complex_add(sec1,sec2);
			tmp3 = gsl_complex_add(tmp3,sec3);

			gsl_matrix_complex_set(h_h,n,k,tmp3);
		}
	}

	x_h_conju = GSLMatrixFree(x_h_conju);
	h_h_conju = GSLMatrixFree(h_h_conju);


}
//------------------culc_eta------------------------------------
void culc_eta()
{
	int n,k,t;
	double tmp1,tmp2;

	for (n = 0; n < N; ++n)
	{
		for (k = 0; k < K; ++k)
		{
			tmp1 = 1;
			tmp2 = 0;
			for (t = 0; t < T; ++t)
			{
				tmp2 += gsl_matrix_get(x_h_abs2,k,t);
				tmp2 = tmp2/(GaussianNoise(sqrt(N0)) + gsl_matrix_get(zeta,n,t));
			}
			tmp1 += tmp2/(double)N;

			gsl_matrix_set(eta,n,k,1/tmp1);
		}
	}
}
//------------------通信路推定器----------------------------------
void channel_estimation()
{
	culc_z();
	culc_zata();
	culc_I();
	
	culc_h_h();
	culc_eta();
	
}
//---------------------------------------------------------------
int main()
{
	int i;
	init();
	// PrintMatrix(stdout,N,T,x_h);
	//x^2の計算
	abs2_matrix(x_h,x_h_abs2,K,T);
	for (i = 0; i < 100; ++i)
	{
		channel_estimation();
		// PrintMatrix(stdout,N,K,h_h);
	}
	PrintMatrix(stdout,N,K,h_h);

	PrintMatrix(stdout,N,K,h);

	// PrintMatrix(stdout,K,T,x_h);

	// PrintMatrix(stdout,K,T,x);

	finish();
	
	return 0;
}