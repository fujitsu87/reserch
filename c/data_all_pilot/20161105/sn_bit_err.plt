set logscale y
set yrange [0.00001:1]
set xlabel "SNR dB"
set ylabel "BIT ERROR RATE"
plot "sn_bit_err_N128_K32_T500_Tp200.dat" w l title "allpilot"
#replot "../data_all_pilot1029/sn_bit_err_N128_K32_T1000_Tp400_allpilot.dat" w l title "allpilot"