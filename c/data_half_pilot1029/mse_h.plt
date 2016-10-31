set logscale y
set yrange [0.001:40]
set xlabel "n"
set ylabel "mse"
plot "mse_h_N128_K32_T1000_Tp400_SN2.dat" w l title "3.01dB"
replot "mse_h_N128_K32_T1000_Tp400_SN4.dat" w l title "6.02dB"
replot "mse_h_N128_K32_T1000_Tp400_SN6.dat" w l
replot "mse_h_N128_K32_T1000_Tp400_SN8.dat" w l title "9.03dB"
replot "mse_h_N128_K32_T1000_Tp400_SN10.dat" w l title "10.00dB"
replot "mse_h_N128_K32_T1000_Tp400_SN20.dat" w l title "13.01dB"
replot "mse_h_N128_K32_T1000_Tp400_SN30.dat" w l 
replot "mse_h_N128_K32_T1000_Tp400_SN40.dat" w l title "16.02dB"
replot "mse_h_N128_K32_T1000_Tp400_SN50.dat" w l
replot "mse_h_N128_K32_T1000_Tp400_SN60.dat" w l title "17.78dB"
