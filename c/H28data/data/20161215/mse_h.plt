set logscale y
set yrange [0.005:1]
set xlabel "n"
set ylabel "MSE"
set key outside
plot "mse_h_N256_K32_T1000_Tp300_SN3.dat" w l title "SNR=0.47dB"
replot "mse_h_N256_K32_T1000_Tp300_SN4.dat" w l title "SNR=0.60dB"
replot "mse_h_N256_K32_T1000_Tp300_SN8.dat" w l title "SNR=0.90dB"
replot "mse_h_N256_K32_T1000_Tp300_SN10.dat" w l title "SNR=1.00dB"
replot "mse_h_N256_K32_T1000_Tp300_SN20.dat" w l title "SNR=1.30dB"
replot "mse_h_N256_K32_T1000_Tp300_SN30.dat" w l title "SNR=1.48dB"
replot "mse_h_N256_K32_T1000_Tp300_SN40.dat" w l title "SNR=1.60dB"
replot "mse_h_N256_K32_T1000_Tp300_SN50.dat" w l title "SNR=1.70dB"
replot "mse_h_N256_K32_T1000_Tp300_SN60.dat" w l title "SNR=1.77dB"
