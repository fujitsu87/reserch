set logscale y
set yrange [0.0001:1]
set xlabel "n"
set ylabel "BIT ERROR RATE"
set key outside
plot "bit_err_N128_K32_T600_Tp250_SN2.dat" w l title "3.01dB"
replot "bit_err_N128_K32_T600_Tp250_SN4.dat" w l title "6.02dB"
replot "bit_err_N128_K32_T600_Tp250_SN6.dat" w l title "7.78dB"
replot "bit_err_N128_K32_T600_Tp250_SN8.dat" w l title "9.03dB"
replot "bit_err_N128_K32_T600_Tp250_SN10.dat" w l title "10.00dB"
replot "bit_err_N128_K32_T600_Tp250_SN20.dat" w l title "13.01dB"
replot "bit_err_N128_K32_T600_Tp250_SN40.dat" w l title "16.02dB"
replot "bit_err_N128_K32_T600_Tp250_SN50.dat" w l title "16.87dB"

