set logscale y
set xlabel "SNR"
set ylabel "MSE"
plot "sn_mes_N256_K32_T128_Tp32.dat" w lp lw 5 title "subspace"
replot "amp.dat" w lp lw 5 title "AMP"
