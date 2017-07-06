set logscale y
set yrange [0.0001:1]
set xlabel "n"
set ylabel "BIT ERROR RATE"
plot "bit_err_N128_K32_T1000_Tp400_SN2.dat" w l title "3.01dB"
replot "bit_err_N128_K32_T1000_Tp400_SN4.dat" w l title "6.02dB"
replot "bit_err_N128_K32_T1000_Tp400_SN8.dat" w l title "9.03dB"
replot "bit_err_N128_K32_T1000_Tp400_SN10.dat" w l title "10.00dB"
replot "bit_err_N128_K32_T1000_Tp400_SN20.dat" w l title "13.01dB"
replot "bit_err_N128_K32_T1000_Tp400_SN40.dat" w l title "16.02dB"
replot "bit_err_N128_K32_T1000_Tp400_SN60.dat" w l title "17.78dB"
