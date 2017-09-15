set logscale y
set xlabel "n"
set ylabel "BER"
plot "mse_h_N256_K32_T64_Tp32_SN10.dat" w l title "p = 0.1"
replot "../20170723_2/mse_h_N256_K32_T64_Tp32_SN10.dat" w l title "p = 0.2"