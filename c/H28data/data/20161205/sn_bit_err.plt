set logscale y
set xlabel "SNR dB"
set ylabel "BER"
set yrange [0.0001:1]
set key outside
set format y "10^{%L}"
plot "sn_bit_err_N256_K32_T1000_Tp300.dat" with linespoints lw 3 notitle
replot "../20161215/sn_bit_err_N256_K32_T1000_Tp300.dat" with linespoints lw 3 notitle