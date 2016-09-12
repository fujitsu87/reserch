#ifndef FFT_C
#define FFT_C 1

#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_vector.h>
#include <math.h>

//----------------複素ベクトルのDFT----------------
void DFT(gsl_vector_complex *x, gsl_vector_complex *y)
{
	gsl_vector_complex_memcpy(y, x);

	if( (y->size & (y->size - 1) ) )	//2の累乗かチェック
	{
		//2の累乗以外
		gsl_fft_complex_wavetable * wavetable = gsl_fft_complex_wavetable_alloc(y->size);
		gsl_fft_complex_workspace * workspace = gsl_fft_complex_workspace_alloc(y->size);

		gsl_fft_complex_forward(y->data, y->stride, y->size, wavetable, workspace);

		gsl_fft_complex_wavetable_free(wavetable);
		gsl_fft_complex_workspace_free(workspace);
	}else{
		//2の累乗
		gsl_fft_complex_radix2_forward(y->data, y->stride, y->size);
	}

	{
		gsl_complex f;
		GSL_SET_COMPLEX(&f, 1.0 / sqrt(y->size), 0);
		gsl_vector_complex_scale(y, f);
	}
}

//----------------複素ベクトルのIDFT----------------
void IDFT(gsl_vector_complex *x, gsl_vector_complex *y)
{
	gsl_vector_complex_memcpy(y, x);

	if( (y->size & (y->size - 1) ) )	//2の累乗かチェック
	{
		//2の累乗以外
		gsl_fft_complex_wavetable * wavetable = gsl_fft_complex_wavetable_alloc(y->size);
		gsl_fft_complex_workspace * workspace = gsl_fft_complex_workspace_alloc(y->size);

		gsl_fft_complex_backward(y->data, y->stride, y->size, wavetable, workspace);

		gsl_fft_complex_wavetable_free(wavetable);
		gsl_fft_complex_workspace_free(workspace);

	}else{
		//2の累乗
		gsl_fft_complex_radix2_backward(y->data, y->stride, y->size);
	}

	{
		gsl_complex f;
		GSL_SET_COMPLEX(&f, 1.0 / sqrt(y->size), 0);
		gsl_vector_complex_scale(y, f);
	}
}

#endif