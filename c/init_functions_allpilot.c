//基地局のアンテナの数
#define N 128
//ユーザ数
#define K 32
//離散時間ステップ数
#define T 600
//pilot信号時間長さ
#define Tp 250
//ユーザ分割数
#define DK 2
//pilot信号の値
const double Pk = 1.0;

//ノイズ　10dB
double N0;
//信号が送られる確率
const double rho_k = 1.0;

double a_k;

int Repeat_flg = 0;

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

gsl_matrix* pilot;



//--------------------------------初期化関数--------------------------------------

void init_h()
{
	int n,k;
	double sigma = 0.5;
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

void init_pilot()
{
	int k,t;
	for (k = 0; k < K; ++k)
	{
		for (t = 0; t < T; ++t)
		{
			
			
			
			// 左のpilot信号
			if(t < Tp)
			{
				gsl_matrix_set(pilot,k,t,1);
			}
			else{
				gsl_matrix_set(pilot,k,t,0);
			}
			
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

void init_xi()
{
	int i,j;
	for (i = 0; i < K; ++i)
	{
		for (j = 0; j < T; ++j)
		{
			double tmp[2] = {1.0,0.0};
			int index = (int)gsl_matrix_get(pilot,i,j);
			gsl_matrix_set(xi ,i,j,tmp[index]);
			gsl_matrix_set(xi_b ,i,j,1);
		}
	}

}

void init_eta()
{
	int i,j;
	for (i = 0; i < N; ++i)
	{
		for (j = 0; j < K; ++j)
		{
			gsl_matrix_set(eta ,i,j,1);
		}
	}
}

void init_x_h()
{
	int k,t;
	gsl_complex def;
	def.dat[0] = 0;
	def.dat[1] = 0;
	for (k = 0; k < K; ++k)
	{
		for (t = 0; t < T; ++t)
		{
			if(gsl_matrix_get(pilot,k,t) == 1)
			{
				gsl_matrix_complex_set(x_h,k,t,gsl_matrix_complex_get(x,k,t));
				gsl_matrix_complex_set(x_b,k,t,gsl_matrix_complex_get(x,k,t));
			}
			else 
			{
				gsl_matrix_complex_set(x_h,k,t,def);
				gsl_matrix_complex_set(x_b,k,t,def);
			}
		}
	}
}

void init(double sn)
{
	int i,j;
	//ノイズ　10dB
	N0 = Pk/sn;
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

	pilot = gsl_matrix_calloc(K,T);

	a_k = sqrt(Pk/(2*rho_k));
	
	//乱数初期化
	RandomNumberInitialization(2);
	
	//通信路h初期化
	init_h();
	
	//pilot信号の初期化　pilotの場合...0 他...1
	init_pilot();

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

	init_x_h();
	// init xi
	init_xi();
	//init eta
	init_eta();

	// PrintMatrix(stdout,K,T,x);
	// PrintMatrix(stdout,K,T,x_h);

}