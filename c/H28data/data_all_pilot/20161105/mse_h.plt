set logscale y
set yrange [0.001:40]
set xlabel "n"
set ylabel "mse"
set key outside
plot "mse_h_N128_K32_T500_Tp200_SN2.dat" w l title "3.01dB"
replot "mse_h_N128_K32_T500_Tp200_SN4.dat" w l title "6.02dB"
replot "mse_h_N128_K32_T500_Tp200_SN6.dat" w l title "7.78dB"
replot "mse_h_N128_K32_T500_Tp200_SN8.dat" w l title "9.03dB"
replot "mse_h_N128_K32_T500_Tp200_SN10.dat" w l title "10.00dB"
replot "mse_h_N128_K32_T500_Tp200_SN20.dat" w l title "13.01dB"
replot "mse_h_N128_K32_T500_Tp200_SN40.dat" w l title "16.02dB"
replot "mse_h_N128_K32_T500_Tp200_SN50.dat" w l title "16.87dB"
replot "mse_h_N128_K32_T500_Tp200_SN70.dat" w l title "17.78dB"
