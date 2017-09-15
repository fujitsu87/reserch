set logscale y
set xlabel "n"
set ylabel "MSE"
plot "mse_h_N256_K32_T1000_Tp256_SN10.dat" w l linewidth 4 title "Tp256"
replot "mse_h_N256_K32_T1000_Tp200_SN10.dat" w l linewidth 4 title "Tp200"
replot "mse_h_N256_K32_T1000_Tp180_SN10.dat" w l linewidth 4 title "Tp180"
replot "mse_h_N256_K32_T128_Tp64_SN10.dat" w l linewidth 4 title "Tp64 p = 0.1"