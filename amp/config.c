#ifndef __CONFIG_C__
#define __CONFIG_C__

//基地局のアンテナの数
#define N 32
//ユーザ数
#define K 40
//離散時間ステップ数
#define T 128
//pilot信号時間長さ(必ずTP>Kにする)
#define Tp 40
//隣接基地局数
#define DK 2
//アンサンブル平均回数
#define ENSEMBLE 1

//反復回数 おおまわり
#define BIG_LOOP 2
//反復回数　小さいループ
#define H_LOOP 10
#define X_LOOP 10

//電力
const double Pk = 1.0;

//pilotの置き方　0...shift 1...sync 2...contamination
const int Pilot_flg = 1;

//ノイズ
double N0;
//信号が送られる確率
const double rho_k = 1.0;
double a_k;

int Repeat_flg = 1;

//ダンピング係数
double a_x= 1.0;
double a_h= 1.0;

//セル間電力差
double cell_diff = 1.0;

gsl_matrix_complex* x;
gsl_matrix_complex* h;
gsl_matrix_complex* w;
gsl_matrix_complex* y;

gsl_matrix_complex* z_d;
gsl_matrix_complex* z_c;

gsl_matrix_complex* x_h;
gsl_matrix_complex* x_h_pre;
gsl_matrix_complex* x_b;
gsl_matrix* x_h_abs2;
gsl_matrix* xi;
gsl_matrix* xi_b;

gsl_matrix_complex* h_h;
gsl_matrix_complex* h_h_pre;
gsl_matrix_complex* h_b;
gsl_matrix* h_h_abs2;
gsl_matrix* eta;
gsl_matrix* eta_b;

gsl_matrix* zeta_c;
gsl_matrix* zeta_d;

gsl_matrix* pilot;

//各ユーザーの電力
double* user_p;

#endif
