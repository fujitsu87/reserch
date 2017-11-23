set logscale y
set xlabel "n"
set ylabel "MSE"
plot "mse_h_N256_K8_T1000_Tp300_SN10.dat" w l linewidth 4 title "whole"
replot "mse_h_my_N256_K8_T1000_Tp300_SN10.dat" w l linewidth 4 title "my"
replot "mse_h_other_N256_K8_T1000_Tp300_SN10.dat" w l linewidth 4 title "other"
