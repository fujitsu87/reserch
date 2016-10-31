set logscale y
set yrange [0.00001:1]
set xlabel "S/N dN"
set ylabel "BIT ERROR RATE"
plot "sn_bit_err_N1024_K4_T1000_Tp400_halfpilot.dat" w l
replot "../data_all_pilot1023/sn_bit_err_N1024_K4_T1000_Tp400_allpilot.dat" w l