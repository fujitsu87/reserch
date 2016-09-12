#ifndef HISTGRAM
#define HISTGRAM 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_histogram.h>

#define EPSILON 1e-10

void Histogram(double min,double max,size_t number,FILE *input_fp,FILE *output_fp)
{
	//データ
	double x;
	//データ数
	int n = 0;
	//ヒストグラム
	gsl_histogram *hist;
	int i;

	//階級値の幅
	double d = (max - min)/(double)number;

	//ヒストグラムの確保
	hist = gsl_histogram_alloc(number);

	//下限 min から上限 max までを等分する階級を設定する。
	gsl_histogram_set_ranges_uniform(hist, min, max);

	//ヒストグラムの全ての階級値を零にする。
	gsl_histogram_reset(hist);
	//ヒストグラムに加えていく
	while( fscanf( input_fp, "%lf",&x) != EOF  )
	{

		if((x > min || fabs(x - min) < EPSILON) && x < max)
		{
			int index = x/d - min/d;
			// 型キャストが失敗したとき
			if (x < hist->range[index] || x > hist->range[index + 1] || fabs(x - hist->range[index + 1]) < EPSILON )
			{
				int upper = number ;
				int lower = 0 ;
				//二分法
				while (upper - lower > 1)
				{
					int mid = (upper + lower) / 2 ; 
					if (x > hist->range[mid] || fabs(x - hist->range[mid]) < EPSILON)
					{
						lower = mid ;
					}
					else
					{
						upper = mid ;
					}
				}
				index = lower;
			}
			hist->bin[index] += 1.0;
			++n;
		}
		else
		{
			printf("%g が見つかりました．\n(min <= x <max)となるように下限値と上限値の値を正しく設定してください．\n",x);
			exit(EXIT_FAILURE);
		}
		
	}

	//正規化する
	for (i = 0;  i < number; i++)
	{
		hist->bin[i] = (double)hist->bin[i]/(double)n;

		//階級値
		double tmp = (hist->range[i]+hist->range[i+1])/2.0;

		//確率密度の推定値
		double density = hist->bin[i] / d;
		//ファイル出力
		fprintf(output_fp, "%g %g %g\n",tmp,hist->bin[i],density);
	}

	gsl_histogram_free(hist);
	return;
}
#endif