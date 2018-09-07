#ifndef __CONFIG_C__
#define __CONFIG_C__

//基地局のアンテナの数
#define N 128
//ユーザ数
#define K 56
//離散時間ステップ数
#define T 16
//pilot信号時間長さ(必ずTP>Kにする)
#define Tp 8
//隣接基地局数
#define DK 7
//アンサンブル平均回数
#define ENSEMBLE 100

//反復回数 おおまわり
#define BIG_LOOP 10
//反復回数　小さいループ
#define H_LOOP 10
#define X_LOOP 10

//電力
const double Pk = 1.0;

//pilotの置き方　0...shift 1...sync 2...contamination
const int Pilot_flg = 2;

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
double cell_diff = 0.1;

gsl_matrix_complex* x;
gsl_matrix_complex* h;
gsl_matrix_complex* w;
gsl_matrix_complex* y;

gsl_vector_complex* Z;

gsl_matrix_complex* x_h;
gsl_matrix_complex* x_h_sub;

gsl_matrix_complex* h_h;

gsl_matrix* eta;

gsl_matrix* pilot;

gsl_matrix_complex *s;
gsl_matrix_complex *y_sub;
gsl_vector_complex *inno;

gsl_matrix_complex *h_sub;
gsl_matrix_complex *h_sub_true;

gsl_matrix_complex* R_x;
gsl_matrix_complex* R_a;

gsl_matrix_complex *x_kro;
gsl_matrix_complex *x_kro_conju;

gsl_matrix_complex *D;
gsl_matrix_complex *D_other;
gsl_matrix_complex *D_large;
gsl_matrix_complex *P;
gsl_matrix_complex *G;
gsl_matrix_complex *G_ans;
gsl_vector_complex* g;

//各ユーザーの電力
double* user_p;

#endif
