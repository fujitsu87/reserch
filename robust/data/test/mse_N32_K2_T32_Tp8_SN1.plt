set logscale y
set xlabel "n"
set ylabel "MSE"
plot "mse_h_N32_K2_T32_Tp8_SN1.dat" w lp linewidth 1 title "whole"
replot "mse_h_my_N32_K2_T32_Tp8_SN1.dat" w lp linewidth 1 title "my"
replot "mse_h_other_N32_K2_T32_Tp8_SN1.dat" w lp linewidth 1 title "other"
