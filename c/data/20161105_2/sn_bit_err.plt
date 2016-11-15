set logscale y
set yrange [0.00001:1]
set xlabel "SNR dB"
set ylabel "BIT ERROR RATE"
plot "sn_bit_err_N128_K32_T600_Tp250.dat" w l title "halfpilot"
replot "../../data_all_pilot/20161105_2/sn_bit_err_N128_K32_T600_Tp250.dat" w l title "allpilot"