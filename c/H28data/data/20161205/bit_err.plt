set logscale y
set yrange [0.000001:1]
set xlabel "n"
set ylabel "BER"
set key outside
set format y "10^{%L}"
plot "bit_err_N256_K32_T1000_Tp300_SN3.dat" w l lw 3 title  "SNR=4.77dB"
replot "bit_err_N256_K32_T1000_Tp300_SN4.dat" w l lw 3 title  "SNR=6.02dB"
replot "bit_err_N256_K32_T1000_Tp300_SN5.dat" w l lw 3 title  "SNR=6.99dB"
replot "bit_err_N256_K32_T1000_Tp300_SN6.dat" w l lw 3 title  "SNR=7.78dB"
replot "bit_err_N256_K32_T1000_Tp300_SN7.dat" w l lw 3 title  "SNR=8.45dB"
replot "bit_err_N256_K32_T1000_Tp300_SN8.dat" w l lw 3 title  "SNR=9.03dB"
replot "bit_err_N256_K32_T1000_Tp300_SN10.dat" w l lw 3 title  "SNR=10.00dB"
replot "bit_err_N256_K32_T1000_Tp300_SN12.dat" w l lw 3 title  "SNR=10.79dB"
replot "bit_err_N256_K32_T1000_Tp300_SN15.dat" w l lw 3 title  "SNR=11.76dB"
replot "bit_err_N256_K32_T1000_Tp300_SN20.dat" w l lw 3 title  "SNR=13.01dB"
