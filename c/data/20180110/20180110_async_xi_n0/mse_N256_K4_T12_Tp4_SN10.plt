set logscale y
set xlabel "n"
set ylabel "MSE"
plot "mse_h_N256_K4_T12_Tp4_SN10.dat" w lp linewidth 1 title "whole"
replot "mse_h_my_N256_K4_T12_Tp4_SN10.dat" w lp linewidth 1 title "my"
replot "mse_h_other_N256_K4_T12_Tp4_SN10.dat" w lp linewidth 1 title "other"
