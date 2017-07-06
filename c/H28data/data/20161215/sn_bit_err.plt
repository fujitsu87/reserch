set logscale y
set xlabel "SNR dB"
set ylabel "BER"
set yrange [0.0001:1]
set key outside
plot "sn_bit_err_N256_K32_T1000_Tp300.dat" with linespoints lw 3 title "no pilot shift"
replot "../20161205/sn_bit_err_N256_K32_T1000_Tp300.dat" with linespoints lw 3 title "pilot shitf"