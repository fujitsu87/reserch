set logscale y
set xlabel "n"
set ylabel "MSE"
plot "mse_h_N256_K16_T64_Tp16_SN10.dat" w l linewidth 4 title "whole"
replot "mse_h_my_N256_K16_T64_Tp16_SN10.dat" w l linewidth 4 title "my user"
replot "mse_h_other_N256_K16_T64_Tp16_SN10.dat" w l linewidth 4 title "other user"