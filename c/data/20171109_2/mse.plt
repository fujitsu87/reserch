set logscale y
set xlabel "n"
set ylabel "MSE"
plot "mse_h_N128_K64_T1024_Tp256_SN10.dat" w l linewidth 4 title "whole"
replot "mse_h_my_N128_K64_T1024_Tp256_SN10.dat" w l linewidth 4 title "my user"
replot "mse_h_other_N128_K64_T1024_Tp256_SN10.dat" w l linewidth 4 title "other user"