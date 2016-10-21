
//AとBの平均2乗誤差を出力
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

double culc_bit_error_rate()
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
				//ビットを間違えた場合
				if(gsl_matrix_complex_get(x_h,k,t).dat[0]*gsl_matrix_complex_get(x,k,t).dat[0] < 0 )
				{
					++count;
				}
				if(gsl_matrix_complex_get(x_h,k,t).dat[1]*gsl_matrix_complex_get(x,k,t).dat[1] < 0 )
				{
					++count;
				}
				num = num + 2;
			}
		}
	}
	ans = (double)count/(double)num;
	return ans;
}