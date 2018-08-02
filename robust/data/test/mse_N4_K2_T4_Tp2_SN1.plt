set logscale y
set xlabel "n"
set ylabel "MSE"
plot "mse_h_N4_K2_T4_Tp2_SN1.dat" w lp linewidth 1 title "whole"
replot "mse_h_my_N4_K2_T4_Tp2_SN1.dat" w lp linewidth 1 title "my"
replot "mse_h_other_N4_K2_T4_Tp2_SN1.dat" w lp linewidth 1 title "other"
