set logscale y
set xlabel "SNR in dB"
set ylabel "BER"
plot "sn_bit_err_N32_K32_T128_Tp64.dat" with points pointtype 7
#replot "../20171116/sn_bit_err_N32_K32_T128_Tp64.dat" with points pointtype 7

