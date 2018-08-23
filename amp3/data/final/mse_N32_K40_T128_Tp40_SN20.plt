set logscale y
set xlabel "n"
set ylabel "MSE"
plot "mse_h_N32_K40_T128_Tp40_SN20.dat" w lp linewidth 1 title "whole"
replot "mse_h_my_N32_K40_T128_Tp40_SN20.dat" w lp linewidth 1 title "my"
replot "mse_h_other_N32_K40_T128_Tp40_SN20.dat" w lp linewidth 1 title "other"
