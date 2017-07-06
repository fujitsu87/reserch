set logscale y
set yrange [0.0001:1]
set xlabel "n"
set ylabel "BER"
set key outside
plot "bit_err_N256_K32_T1000_Tp300_SN3.dat" w l title "SNR=0.47dB"
replot "bit_err_N256_K32_T1000_Tp300_SN4.dat" w l title "SNR=0.60dB"
replot "bit_err_N256_K32_T1000_Tp300_SN8.dat" w l title "SNR=0.90dB"
replot "bit_err_N256_K32_T1000_Tp300_SN10.dat" w l title "SNR=1.00dB"
replot "bit_err_N256_K32_T1000_Tp300_SN20.dat" w l title "SNR=1.30dB"
#replot "bit_err_N256_K32_T1000_Tp300_SN30.dat" w l title "SNR=1.48dB"
#replot "bit_err_N256_K32_T1000_Tp300_SN40.dat" w l title "SNR=1.60dB"
#replot "bit_err_N256_K32_T1000_Tp300_SN50.dat" w l title "SNR=1.70dB"
#replot "bit_err_N256_K32_T1000_Tp300_SN60.dat" w l title "SNR=1.77dB"
