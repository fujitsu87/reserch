set logscale y
set xlabel "n"
set ylabel "MSE"
plot "mse_h_N16_K32_T128_Tp32_SN10.dat" w lp linewidth 1 title "whole"
replot "mse_h_my_N16_K32_T128_Tp32_SN10.dat" w lp linewidth 1 title "my"
replot "mse_h_other_N16_K32_T128_Tp32_SN10.dat" w lp linewidth 1 title "other"