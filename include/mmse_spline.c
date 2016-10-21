#ifndef BPSKMMSESPLINE
#define  BPSKMMSESPLINE 1


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "../include/mmse.c"

//----------global variable--------------------
#define NDATA 185
double BPSKMMSESpline_data[NDATA] =
{
	-80,
	-50,
	-30,
	-20,
	-15,
	-14,
	-13,
	-12,
	-11,
	-10,
	-9,
	-8.5,
	-8,
	-7.5,
	-7,
	-6.5,
	-6,
	-5.5,
	-5,
	-4.5,
	-4,
	-3.6,
	-3.5,
	-3.4,
	-3.3,
	-3.2,
	-3.1,
	-3,
	-2.95,
	-2.9,
	-2.88,
	-2.87,
	-2.85,
	-2.8,
	-2.75,
	-2.5,
	-2.25,
	-2,
	-1.75,
	-1.5,
	-1.25,
	-1,
	-0.75,
	-0.5,
	-0.25,
	0,
	0.25,
	0.5,
	0.6,
	0.7,
	0.75,
	1,
	1.25,
	1.5,
	1.75,
	2,
	2.25,
	2.5,
	2.6,
	2.7,
	2.8,
	2.9,
	3,
	3.1,
	3.2,
	3.3,
	3.4,
	3.5,
	3.6,
	3.7,
	3.8,
	3.9,
	4,
	4.1,
	4.2,
	4.3,
	4.4,
	4.5,
	4.615,
	4.75,
	4.815,
	5,
	5.125,
	5.25,
	5.315,
	5.5,
	5.615,
	5.75,
	5.825,
	6,
	6.125,
	6.25,
	6.375,
	6.4,
	6.5,
	6.615,
	6.75,
	6.9,
	7,
	7.125,
	7.25,
	7.3,
	7.4,
	7.5,
	7.6,
	7.75,
	7.8,
	7.9,
	8,
	8.1,
	8.2,
	8.3,
	8.4,
	8.5,
	8.6,
	8.7,
	8.8,
	8.9,
	9,
	9.5,
	9.6,
	9.7,
	9.8,
	9.9,
	10,
	10.25,
	10.5,
	10.75,
	11,
	11.25,
	11.5,
	11.75,
	12,
	12.5,
	12.75,
	13,
	13.5,
	14,
	14.5,
	15,
	16,
	17,
	18,
	19,
	20,
	21,
	22,
	24,
	25,
	26,
	27,
	28,
	29,
	30,
	35,
	36,
	37,
	38,
	39,
	40,
	41,
	42,
	43,
	44,
	45,
	46,
	47,
	48,
	49,
	50,
	60,
	65,
	66,
	67,
	68,
	69,
	70,
	72,
	73,
	74,
	75,
	77,
	78,
	79,
	80
};

gsl_interp_accel *BPSKMMSESpline_acc;
gsl_spline *BPSKMMSESpline_spline;

//-----------function to be complemented---------------
double interpolation_func(double x)
{
	return BPSKMMSE(x);
}
//---------initialization of data----------------------
void init_data(double *x,double *y)
{
	int i;
	for(i = 0; i < NDATA; ++i)
	{
		x[i] = pow(10,BPSKMMSESpline_data[i]/10);
		y[i] = interpolation_func(x[i]);
	}
	return;
}
//----------initialization of the spline function-------------
void SplineInit()
{
	double x[NDATA];
	double y[NDATA];
	init_data(x,y);

	BPSKMMSESpline_acc = gsl_interp_accel_alloc();
	BPSKMMSESpline_spline = gsl_spline_alloc(gsl_interp_cspline,NDATA);
	gsl_spline_init(BPSKMMSESpline_spline, x, y, NDATA);

}
//-------------return complemented value-------------
double  BPSKMMSESpline(double snr)
{
	return  gsl_spline_eval(BPSKMMSESpline_spline, snr, BPSKMMSESpline_acc);
}
#endif