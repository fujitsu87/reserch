set logscale y
unset key
set yrange [0.001:40]
set xlabel "n"
set ylabel "mse"
plot "mse_h_N1024_K32_T1000_Tp400_SN2_half.dat" w l
replot "mse_h_N1024_K32_T1000_Tp400_SN4_half.dat" w l
replot "mse_h_N1024_K32_T1000_Tp400_SN6_half.dat" w l
replot "mse_h_N1024_K32_T1000_Tp400_SN8_half.dat" w l
replot "mse_h_N1024_K32_T1000_Tp400_SN10_half.dat" w l
replot "mse_h_N1024_K32_T1000_Tp400_SN20_half.dat" w l
replot "mse_h_N1024_K32_T1000_Tp400_SN30_half.dat" w l
replot "mse_h_N1024_K32_T1000_Tp400_SN40_half.dat" w l
replot "mse_h_N1024_K32_T1000_Tp400_SN50_half.dat" w l
replot "mse_h_N1024_K32_T1000_Tp400_SN60_half.dat" w l
replot "mse_h_N1024_K32_T1000_Tp400_SN70_half.dat" w l
replot "mse_h_N1024_K32_T1000_Tp400_SN80_half.dat" w l
replot "mse_h_N1024_K32_T1000_Tp400_SN90_half.dat" w l
replot "mse_h_N1024_K32_T1000_Tp400_SN100_half.dat" w l