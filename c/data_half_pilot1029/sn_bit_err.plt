set logscale y
set yrange [0.00001:1]
set xlabel "S/N dB"
set ylabel "BIT ERROR RATE"
plot "sn_bit_err_N128_K32_T1000_Tp400_halfpilot.dat" w l title "halfpilot"
replot "../data_all_pilot1029/sn_bit_err_N128_K32_T1000_Tp400_allpilot.dat" w l title "allpilot"