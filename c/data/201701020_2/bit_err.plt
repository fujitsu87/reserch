set logscale y
set xlabel "n"
set ylabel "MSE"
plot "mse_h_N256_K48_T384_Tp128_SN10.dat" w l linewidth 4 title "a=0.1 p=0.1"
#replot "./2/mse_h_N256_K32_T128_Tp64_SN10.dat" w l linewidth 4 title "a=0.2 p=1"
#replot "./1/mse_h_N256_K32_T128_Tp64_SN10.dat" w l linewidth 4 title "Tp64 p=0.1"