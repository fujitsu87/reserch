set logscale y
set yrange [0.00001:1]
set xlabel "S/N dB"
set ylabel "BIT ERROR RATE"
plot "sn_bit_err_N2000K8T1000Tp400R3_100_100_halfPilot.txt" w l
replot "../data_all_pilot/sn_bit_err_N2000K8T1000Tp400R3_100_100_allPilot.txt" w l