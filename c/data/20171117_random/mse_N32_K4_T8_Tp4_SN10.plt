set logscale y
set xlabel "n"
set ylabel "MSE"
plot "mse_h_N32_K4_T8_Tp4_SN10.dat" w lp linewidth 4 title "whole"
replot "mse_h_my_N32_K4_T8_Tp4_SN10.dat" w lp linewidth 4 title "my"
replot "mse_h_other_N32_K4_T8_Tp4_SN10.dat" w lp linewidth 4 title "other"
