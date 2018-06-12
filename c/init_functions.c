#include "./config.c"

//--------------------------------初期化関数--------------------------------------

void init_h()
{
	int n,k;
	double sigma = 1;
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
			if(Pilot_flg == 1){
				// 左のpilot信号
				if(t < Tp)
				{
					gsl_matrix_set(pilot,k,t,1);
				}
				else{
					gsl_matrix_set(pilot,k,t,0);
				}
			}
			else if(Pilot_flg == 2){
				// 左のpilot信号
				if(t < Tp)
				{
					gsl_matrix_set(pilot,k,t,1);
				}
				else{
					gsl_matrix_set(pilot,k,t,0);
				}
			}
			else {
				//左下のpilot信号
				if(k >= K/DK && t < Tp)
				{
					gsl_matrix_set(pilot,k,t,1);
				}
				//右上のpilot信号
				else if(k < K/DK && t >= T-Tp)
				{
					gsl_matrix_set(pilot,k,t,1);
				}
				else{
					gsl_matrix_set(pilot,k,t,0);
				}
/* 基地局が3つの場合
				//左下のpilot信号
				if(k >= 2*K/DK && t < Tp)
				{
					gsl_matrix_set(pilot,k,t,1);
				}
				//下から二段目のpilot信号
				else if(k >= K/DK && k < 2*K/DK && t >= Tp && t < 2*Tp)
				{
					gsl_matrix_set(pilot,k,t,1);
				}
				//右上のpilot信号
				else if(k < K/DK && t >= T-Tp)
				{
					gsl_matrix_set(pilot,k,t,1);
				}
				else{
					gsl_matrix_set(pilot,k,t,0);
				}
*/
			}
			
		}
	}
	// PrintRealMatrix(stdout,K,T,pilot);
}

int decision_hadamard_size()
{
	int n = 0;
	int matrix_size = 1;
	while(1)
	{
		if(matrix_size >= Tp)
		{
			return n;
		}
		matrix_size = matrix_size*2;
		++n;
	}
}

void init_x()
{
	int k,t,i;
	int n = decision_hadamard_size();
	int size = (int)pow(2,n);
	int *P = (int *)malloc(sizeof(int)*size);

	gsl_matrix* orthogonal = gsl_matrix_calloc(size,size);
	gsl_matrix* re_pilot = gsl_matrix_calloc(size,size);
	gsl_matrix* im_pilot = gsl_matrix_calloc(size,size);
	orthogonal = hadamard(n);
	
	Permutation(size,P);
	GslMatrixRealInterleaver(size,size,orthogonal,re_pilot,P);
	Permutation(size,P);
	GslMatrixRealInterleaver(size,size,orthogonal,im_pilot,P);
	free(orthogonal);
	free(P);

	//for data 
	int data_size = T-Tp;
	double *data = (double *)malloc(sizeof(double)*data_size);
	double *data_tmp = (double *)malloc(sizeof(double)*data_size);
	for(i=0;i<data_size;i++)
	{
		if(i < data_size*rho_k) data[i] = 1.0;
		else data[i] = 0.0;
	}
	P = (int *)malloc(sizeof(int)*data_size);

	for (k = 0; k < K; ++k)
	{
		double p_p = sqrt(user_p[k]/2);
		a_k = sqrt(user_p[k]/(2*rho_k));

		Permutation(data_size,P);
		RealInterleaver(data_size,data,data_tmp,P);
		for (t = 0; t < T; ++t)
		{
			double tmp[2] = {-1.0,1.0};
			gsl_complex z;

			//左のpilot 同期パイロットの場合
			if(Pilot_flg == 1 && t < Tp && gsl_matrix_get(pilot,k,t)==1)
			{
				// z = gsl_matrix_complex_get(x,k,T-1-t);
				// 左のpilot信号
				if(t < Tp)
				{
					// GSL_SET_COMPLEX(&z,a_k*tmp[UniformBit()],a_k*tmp[UniformBit()]);
					GSL_SET_COMPLEX(&z,p_p*gsl_matrix_get(re_pilot,k,t),p_p*gsl_matrix_get(im_pilot,k,t));
				}
	
			}
/*			
			//右上のpilot信号
			else if(Pilot_flg == 0 && k < K/DK && t >= T-Tp)
			{
				GSL_SET_COMPLEX(&z,p_p*gsl_matrix_get(re_pilot,k,t-(T-Tp)),p_p*gsl_matrix_get(im_pilot,k,t-(T-Tp)));
			}
			//左下
			else if(Pilot_flg == 0 && k >= 2*K/DK && t < Tp)
			{
				GSL_SET_COMPLEX(&z,sqrt(cell_diff)*p_p*gsl_matrix_get(re_pilot,k-K/DK,t),sqrt(cell_diff)*p_p*gsl_matrix_get(im_pilot,k-K/DK,t));
			}
			//下から二段目のpilot信号
			else if(Pilot_flg == 0 && k >= K/DK && k < 2*K/DK && t >= Tp && t < 2*Tp)
			{
				GSL_SET_COMPLEX(&z,sqrt(cell_diff)*p_p*gsl_matrix_get(re_pilot,k-K/DK,t-Tp),sqrt(cell_diff)*p_p*gsl_matrix_get(im_pilot,k-K/DK,t-Tp));
			}
*/
// DK=2のとき
			//右上のpilot信号
			else if(Pilot_flg == 0 && k < K/DK && gsl_matrix_get(pilot,k,t)==1)
			{
				//アダマール
				// GSL_SET_COMPLEX(&z,p_p*gsl_matrix_get(re_pilot,k,t-(T-Tp)),p_p*gsl_matrix_get(im_pilot,k,t-(T-Tp)));
				
				//乱数
				GSL_SET_COMPLEX(&z,p_p*tmp[UniformBit()],p_p*tmp[UniformBit()]);
			}

			//左下のpilot信号　シフトパイロットの場合
			else if(Pilot_flg == 0 && k >= K/DK && gsl_matrix_get(pilot,k,t)==1)
			{
				//乱数
				GSL_SET_COMPLEX(&z,p_p*gsl_matrix_get(re_pilot,k-K/DK,t),p_p*gsl_matrix_get(im_pilot,k-K/DK,t));
				
				// z = gsl_matrix_complex_get(x,k-K/DK,t+T-Tp);
				z.dat[0] = sqrt(cell_diff) * z.dat[0];
				z.dat[1] = sqrt(cell_diff) * z.dat[1];
			}

			//下のpilot信号 contaminationの場合
			else if(Pilot_flg == 2 && k >= K/DK)
			{
				z = gsl_matrix_complex_get(x,k-K/DK,t);
			}
			else 
			{
/*
				if(rho_k > gsl_rng_uniform(RAN)){
					if(k < K/DK){
						GSL_SET_COMPLEX(&z,a_k*tmp[UniformBit()],a_k*tmp[UniformBit()]);
					}
					else{
						GSL_SET_COMPLEX(&z,sqrt(cell_diff)*a_k*tmp[UniformBit()],sqrt(cell_diff)*a_k*tmp[UniformBit()]);
					}
				}
				else GSL_SET_COMPLEX(&z,0,0);
*/
				//permutationする
				//sync
				if(Pilot_flg == 1)
				{
						GSL_SET_COMPLEX(&z,
							a_k*data_tmp[t-Tp]*tmp[UniformBit()],
							a_k*data_tmp[t-Tp]*tmp[UniformBit()]);
				}
				//shift
				else{
					if(k < K/DK){
						GSL_SET_COMPLEX(&z,
							a_k*data_tmp[t]*tmp[UniformBit()],
							a_k*data_tmp[t]*tmp[UniformBit()]);
					}
					else{
						GSL_SET_COMPLEX(&z,
							a_k*data_tmp[t-Tp]*tmp[UniformBit()],
							a_k*data_tmp[t-Tp]*tmp[UniformBit()]);
					}
				}
			}

			gsl_matrix_complex_set(x,k,t,z);
		}
	}
	free(re_pilot);
	free(im_pilot);
	free(P);
	free(data);
	free(data_tmp);

	// PrintRealMatrix(stdout,K,T,pilot);
	// PrintRealMatrix(stdout,size,size,im_pilot);
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

void init_user_p()
{
	int k;
	for (k = 0; k < K; ++k)
	{
		//自セルユーザーの場合
		if(k < K/DK){
			user_p[k] = Pk + GaussianNoise(0);
			// user_p[k] = Pk + GaussianNoise(N0*0.1);
		}
		//他セルユーザーの場合
		else{
			user_p[k] = cell_diff*Pk + GaussianNoise(cell_diff*0);
			// user_p[k] = cell_diff*Pk + GaussianNoise(cell_diff*N0*0.1);
		}

		if(user_p[k] < 0)user_p[k]=0.0;
	}

}

void init_xi_tmp()
{
	int i,j;
	for (i = 0; i < K; ++i)
	{
		for (j = 0; j < T; ++j)
		{
			double tmp[2] = {Pk,0.0};
			int index = (int)gsl_matrix_get(pilot,i,j);
			// if(j >= Tp && j < (T-Tp) && Pilot_flg==0)index = 1;
			// else if(Pilot_flg==1)index = 1;
			gsl_matrix_set(xi ,i,j,0.5);
			// gsl_matrix_set(xi ,i,j,Pk);
			gsl_matrix_set(xi_b ,i,j,Pk);
		}
	}
	// PrintRealMatrix(stdout,K,T,xi);
}

void init_xi_first()
{
	int i,j;
	for (i = 0; i < K; ++i)
	{
		for (j = 0; j < T; ++j)
		{
			double tmp[2] = {Pk,0.0};
			int index = (int)gsl_matrix_get(pilot,i,j);
			if(j >= Tp && j < (T-Tp) && Pilot_flg==0)index = 1;
			else if(Pilot_flg==1)index = 1;
			gsl_matrix_set(xi ,i,j,tmp[index]);
			// gsl_matrix_set(xi ,i,j,Pk);
			gsl_matrix_set(xi_b ,i,j,Pk);
			// gsl_matrix_set(xi ,i,j,tmp[index]);
		}
	}
	// PrintRealMatrix(stdout,K,T,xi);
}
void init_xi()
{
	int i,j;
	for (i = 0; i < K; ++i)
	{
		for (j = 0; j < T; ++j)
		{
			double tmp[2] = {Pk,0.0};
			int index = (int)gsl_matrix_get(pilot,i,j);
			// if(j >= Tp && j < (T-Tp))index = 1;
			gsl_matrix_set(xi ,i,j,tmp[index]);
			// gsl_matrix_set(xi ,i,j,Pk);
			gsl_matrix_set(xi_b ,i,j,Pk);
		}
	}
	// PrintRealMatrix(stdout,K,T,xi);
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
	def.dat[0] = 0.0;
	def.dat[1] = 0.0;
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

void init_x_h_tmp()
{
	int k,t,n;
	gsl_complex def,z;
	def.dat[0] = 0;
	def.dat[1] = 0;

	gsl_matrix_complex *HH = gsl_matrix_complex_calloc(K,N);
	gsl_matrix_complex_transpose_memcpy(HH,h);
	for(k=0;k<K;k++){
		for(n=0;n<N;n++){
			if(gsl_matrix_get(pilot,k,t) == 1)
			{
				z=gsl_matrix_complex_get(HH,k,n);
				z=gsl_complex_conjugate(z);
				z=gsl_complex_mul_real(z,1/sqrt(N));
				gsl_matrix_complex_set(HH,k,n,z);
			}
		}
	}

	AB(HH,y,x_b);

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
				// gsl_matrix_complex_set(x_b,k,t,def);
			}
		}
	}
	HH=GSLMatrixFree(HH);
}

void init(double sn)
{
	int i,j;
	//ノイズ　10dB
	
	N0 = Pk/(double)sn;
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

	user_p  = (double *)malloc( sizeof(double) *K );

	a_k = sqrt(Pk/(2*rho_k));

	//各ユーザの電力初期化
	init_user_p();

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

	// PrintMatrix(stdout,N,K,h);
	// PrintMatrix(stdout,K,T,x);
	// PrintMatrix(stdout,N,T,w);
	// PrintMatrix(stdout,N,T,y);

	init_x_h();
	// init xi
	init_xi();
	//init eta
	init_eta();

	// PrintMatrix(stdout,K,T,x);
	// PrintMatrix(stdout,K,T,x_h);

}