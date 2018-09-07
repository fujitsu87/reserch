#ifndef __CONFIG_C__
#define __CONFIG_C__

//基地局のアンテナの数
#define N 128
//ユーザ数
#define K 32
//離散時間ステップ数
#define T 128
//pilot信号時間長さ(必ずTP>Kにする)
#define Tp 32
//隣接基地局数
#define DK 2
//アンサンブル平均回数
#define ENSEMBLE 1

//反復回数 おおまわり
#define BIG_LOOP 10
//反復回数　小さいループ
#define H_LOOP 10
#define X_LOOP 10

//電力
const double Pk = 1.0;

//pilotの置き方　0...shift 1...sync 2...contamination
const int Pilot_flg = 0;

//ノイズ
double N0;
//信号が送られる確率
const double rho_k = 1.0;
double a_k;

int Repeat_flg = 1;

//ダンピング係数
double a_x= 0.5;
double a_h= 0.5;

//セル間電力差
double cell_diff = 1.0;

gsl_matrix_complex* x;
gsl_matrix_complex* h;
gsl_matrix_complex* w;
gsl_matrix_complex* y;

gsl_matrix_complex* z;

gsl_matrix_complex* x_h;

gsl_matrix_complex* h_h;

gsl_matrix* eta;

gsl_matrix* pilot;

gsl_matrix_complex *s;
gsl_matrix_complex *y_sub;

gsl_matrix_complex *h_sub;
gsl_matrix_complex *h_sub_true;

//各ユーザーの電力
double* user_p;

#endif
