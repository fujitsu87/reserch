unset logscale
set yrange [0:1]
set xlabel "n"
set ylabel "BIT ERROR RATE"
unset key
plot "bit_err_N1024_K32_T1000_Tp400_SN2_half.dat" w l
replot "bit_err_N1024_K32_T1000_Tp400_SN4_half.dat" w l
replot "bit_err_N1024_K32_T1000_Tp400_SN6_half.dat" w l
replot "bit_err_N1024_K32_T1000_Tp400_SN8_half.dat" w l
replot "bit_err_N1024_K32_T1000_Tp400_SN10_half.dat" w l
replot "bit_err_N1024_K32_T1000_Tp400_SN20_half.dat" w l
replot "bit_err_N1024_K32_T1000_Tp400_SN30_half.dat" w l
replot "bit_err_N1024_K32_T1000_Tp400_SN40_half.dat" w l
replot "bit_err_N1024_K32_T1000_Tp400_SN50_half.dat" w l
replot "bit_err_N1024_K32_T1000_Tp400_SN60_half.dat" w l
replot "bit_err_N1024_K32_T1000_Tp400_SN70_half.dat" w l
replot "bit_err_N1024_K32_T1000_Tp400_SN80_half.dat" w l
replot "bit_err_N1024_K32_T1000_Tp400_SN90_half.dat" w l
replot "bit_err_N1024_K32_T1000_Tp400_SN100_half.dat" w l