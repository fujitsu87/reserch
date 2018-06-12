
//AとBの平均2乗誤差(MSE)を出力
double culc_complex_mse(int n,int m,gsl_matrix_complex * A,gsl_matrix_complex * B)
{
	int num = n*m;
	int i,j;
	double ans = 0;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
		{
			ans += gsl_complex_abs2(gsl_complex_sub(gsl_matrix_complex_get(A,i,j),gsl_matrix_complex_get(B,i,j)));
		}
	}
	ans = ans/num;
	return ans;
}

//etaと(h_h-h)^2の平均2乗誤差(MSE)を出力
double culc_mse_eta_h()
{
	int num = N*K;
	int i,j;
	double ans = 0,tmp;
	double h_mse,eta_val;
	for (i = 0; i < N; ++i)
	{
		for (j = 0; j < K; ++j)
		{
			h_mse = gsl_complex_abs2(gsl_complex_sub(gsl_matrix_complex_get(h,i,j),gsl_matrix_complex_get(h_h,i,j)));
			eta_val = gsl_matrix_get(eta,i,j);
			tmp = h_mse - eta_val;
			tmp = tmp*tmp;
			ans += tmp;
		}
	}
	ans = ans/num;
	return ans;
}

//Aの平均を計算
double culc_avr(int n,int m,gsl_matrix * A)
{
	int num = n*m;
	int i,j;
	double ans = 0;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
		{
			ans += gsl_matrix_get(A,i,j);
		}
	}
	ans = ans/num;
	return ans;
}

//AとBの平均2乗誤差(MSE)を出力 自セル
double culc_complex_mse_etc(int n,int m,gsl_matrix_complex * A,gsl_matrix_complex * B,int jstart)
{
	int num = n*(m-jstart);
	int i,j;
	double ans = 0;
	for (i = 0; i < n; ++i)
	{
		for (j = jstart; j < m; ++j)
		{
			ans += gsl_complex_abs2(gsl_complex_sub(gsl_matrix_complex_get(A,i,j),gsl_matrix_complex_get(B,i,j)));
		}
	}
	ans = ans/num;
	return ans;
}

//AとBの平均2乗誤差を出力
double culc_complex_mse_matrix(int n,int m,gsl_matrix_complex * A,gsl_matrix_complex * B)
{
	int num = n*m;
	int i,j;
	double ans = 0;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
		{
			printf(" %g ", gsl_complex_abs2(gsl_complex_sub(gsl_matrix_complex_get(A,i,j),gsl_matrix_complex_get(B,i,j))));
		}
		printf("\n");
	}
	ans = ans/num;
	return ans;
}

double culc_bit_error_rate(FILE* fp_bit_err,int kstart,int kend)
{
	int k,t;
	double ans = 0;

	int count = 0;
	int num = 0;
	for (k = kstart; k < kend; ++k)
	{
		for (t = 0; t < T; ++t)
		{	
			//pilot信号でない場合
			if(gsl_matrix_get(pilot,k,t) == 0)
			{
				double re = gsl_matrix_complex_get(x_h,k,t).dat[0];
				double im = gsl_matrix_complex_get(x_h,k,t).dat[1];

				//ゼロ符号の場合 y < -x + a_k
				if( (fabs(re)+fabs(im)) < a_k)
				{
					//本当はゼロではないとき
					if (gsl_matrix_complex_get(x,k,t).dat[0] != 0.0)
					{
						count = count +2;
					}
				}
				else{
					//本当はゼロ符号のとき
					if(gsl_matrix_complex_get(x,k,t).dat[0] == 0.0)
					{
						count = count +2;
					}
					else
					{ 
						if(re*gsl_matrix_complex_get(x,k,t).dat[0] < 0 )
						{
							++count;
						}
						if(im*gsl_matrix_complex_get(x,k,t).dat[1] < 0 )
						{
							++count;
						}
						if(isnan(im*gsl_matrix_complex_get(x,k,t).dat[1])){
							fprintf(fp_bit_err,"exist nan\n");
							exit(0);
						}
					}
				}
				num = num + 2;
			}
		}
	}
	// PrintMatrix(stdout,K,T,x_h);
	// PrintMatrix(stdout,K,T,x);
	ans = (double)count/(double)num;
	return ans;
}

/*
double culc_bit_error_rate(FILE* fp_bit_err)
{
	int k,t;
	double ans = 0;

	int count = 0;
	int num = 0;
	for (k = 0; k < K; ++k)
	{
		for (t = 0; t < T; ++t)
		{	
			//pilot信号でない場合
			if(gsl_matrix_get(pilot,k,t) == 0)
			{
				double re = gsl_matrix_complex_get(x_h,k,t).dat[0];
				double im = gsl_matrix_complex_get(x_h,k,t).dat[1];
				//本当はゼロ符号のとき
				if(gsl_matrix_complex_get(x,k,t).dat[0] == 0.0)
				{
					count = count +2;
				}
				else
				{ 
					if(re*gsl_matrix_complex_get(x,k,t).dat[0] < 0 )
					{
						++count;
					}
					if(im*gsl_matrix_complex_get(x,k,t).dat[1] < 0 )
					{
						++count;
					}
					if(isnan(im*gsl_matrix_complex_get(x,k,t).dat[1])){
						fprintf(fp_bit_err,"exist nan\n");
						exit(0);
					}
				}
				num = num + 2;
			}
		}
	}
	ans = (double)count/(double)num;
	return ans;
}
*/
double culc_abs2_all_element(int n,int m,gsl_matrix_complex * A)
{
	int num = n*m;
	int i,j;
	double ans = 0;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
		{
			ans += gsl_complex_abs2(gsl_matrix_complex_get(A,i,j));
		}
	}
	ans = ans/num;
	return ans;
}

double culc_abs2_all_element_real(int n,int m,gsl_matrix * A)
{
	int num = n*m;
	int i,j;
	double ans = 0;
	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
		{
			ans += gsl_matrix_get(A,i,j)*gsl_matrix_get(A,i,j);
		}
	}
	ans = ans/num;
	return ans;
}

void output_estimatedata(FILE *fp)
{
	int k,t;
	for (k = 0; k < K; ++k)
	{
		for (t = 0; t < T; ++t)
		{
			if(gsl_matrix_get(pilot,k,t) == 0)
				fprintf(fp,"%f %f\n",gsl_matrix_complex_get(x_h,k,t).dat[0],gsl_matrix_complex_get(x_h,k,t).dat[1]);
		}
	}

}

void culc_standard_deviation(double *standard_deviation,double *mse_h_n,double * mse_h_en,double *bit_err_n,double *bit_err_en)
{
	int i;
	 //MSEとBER分散と標準偏差を求める
    double variance[2]={0,0};
    for (i = 0; i < ENSEMBLE; ++i)
    {
        double tmp[2];
        tmp[0] = mse_h_n[BIG_LOOP*H_LOOP-1] - mse_h_en[i];
        tmp[1] = bit_err_n[BIG_LOOP*H_LOOP-1] - bit_err_en[i];

        variance[0] += tmp[0] * tmp[0];
        variance[1] += tmp[1] * tmp[1];
    }
    variance[0] = variance[0] /((double)ENSEMBLE -1);
    variance[1] = variance[1] /((double)ENSEMBLE -1); 
    standard_deviation[0] = sqrt(variance[0]);
    standard_deviation[1] = sqrt(variance[1]);
    // PrintMatrix(stdout,K,T,x_h);

}

void create_plotfile(char *outputpath)
{
	FILE *fp_mse_plt,*fp_ber_plt;
	char mse_h_plt[100];
	char ber_plt[100];
	sprintf(ber_plt,"%sber_N%d_K%d_T%d_Tp%d_SN%.f.plt",outputpath,N,K,T,Tp,Pk/N0);
	sprintf(mse_h_plt,"%smse_N%d_K%d_T%d_Tp%d_SN%.f.plt",outputpath,N,K,T,Tp,Pk/N0);
	
	if ((fp_mse_plt = fopen(mse_h_plt, "w")) == NULL) {
        printf("can not open%s\n",mse_h_plt);
        return;
    }
	if ((fp_ber_plt = fopen(ber_plt, "w")) == NULL) {
        printf("can not open%s\n",ber_plt);
        return;
    }

	fprintf(fp_ber_plt,"set logscale y\n");
	fprintf(fp_ber_plt,"set xlabel \"n\"\n");
	fprintf(fp_ber_plt,"set ylabel \"BER\"\n");
	fprintf(fp_ber_plt,"plot \"bit_err_N%d_K%d_T%d_Tp%d_SN%.f.dat\" w lp linewidth 1 title \"whole\"\n",N,K,T,Tp,Pk/N0);
	fprintf(fp_ber_plt,"replot \"bit_err_my_N%d_K%d_T%d_Tp%d_SN%.f.dat\" w lp linewidth 1 title \"my\"\n",N,K,T,Tp,Pk/N0);
	fprintf(fp_ber_plt,"replot \"bit_err_other_N%d_K%d_T%d_Tp%d_SN%.f.dat\" w lp linewidth 1 title \"other\"\n",N,K,T,Tp,Pk/N0);

	fprintf(fp_mse_plt,"set logscale y\n");
	fprintf(fp_mse_plt,"set xlabel \"n\"\n");
	fprintf(fp_mse_plt,"set ylabel \"MSE\"\n");
	fprintf(fp_mse_plt,"plot \"mse_h_N%d_K%d_T%d_Tp%d_SN%.f.dat\" w lp linewidth 1 title \"whole\"\n",N,K,T,Tp,Pk/N0);
	fprintf(fp_mse_plt,"replot \"mse_h_my_N%d_K%d_T%d_Tp%d_SN%.f.dat\" w lp linewidth 1 title \"my\"\n",N,K,T,Tp,Pk/N0);
	fprintf(fp_mse_plt,"replot \"mse_h_other_N%d_K%d_T%d_Tp%d_SN%.f.dat\" w lp linewidth 1 title \"other\"\n",N,K,T,Tp,Pk/N0);

	fclose(fp_mse_plt);
	fclose(fp_ber_plt);
}