set logscale y
set xlabel "n"
set ylabel "MSE"
plot "mse_h_N125_K32_T128_Tp64_SN10.dat" w lp linewidth 4 title "whole"
replot "mse_h_my_N125_K32_T128_Tp64_SN10.dat" w lp linewidth 4 title "my"
replot "mse_h_other_N125_K32_T128_Tp64_SN10.dat" w lp linewidth 4 title "other"
