set logscale y
set xlabel "n"
set ylabel "MSE"
plot "mse_h_N16_K2_T16_Tp8_SN100.dat" w lp linewidth 1 title "whole"
replot "mse_h_my_N16_K2_T16_Tp8_SN100.dat" w lp linewidth 1 title "my"
replot "mse_h_other_N16_K2_T16_Tp8_SN100.dat" w lp linewidth 1 title "other"
