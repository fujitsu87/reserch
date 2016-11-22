set logscale y
set yrange [0.0001:1]
set xlabel "n"
set ylabel "BIT ERROR RATE"
set key outside
plot "bit_err_N160_K100_T100_Tp40_SN1.dat" w l 
replot "bit_err_N160_K100_T100_Tp40_SN2.dat" w l 
replot "bit_err_N160_K100_T100_Tp40_SN3.dat" w l 
replot "bit_err_N160_K100_T100_Tp40_SN4.dat" w l 
replot "bit_err_N160_K100_T100_Tp40_SN5.dat" w l 
replot "bit_err_N160_K100_T100_Tp40_SN8.dat" w l 
replot "bit_err_N160_K100_T100_Tp40_SN10.dat" w l 
replot "bit_err_N160_K100_T100_Tp40_SN20.dat" w l 
replot "bit_err_N160_K100_T100_Tp40_SN30.dat" w l 
replot "bit_err_N160_K100_T100_Tp40_SN40.dat" w l 
replot "bit_err_N160_K100_T100_Tp40_SN50.dat" w l 
replot "bit_err_N160_K100_T100_Tp40_SN60.dat" w l 
replot "bit_err_N160_K100_T100_Tp40_SN70.dat" w l 
replot "bit_err_N160_K100_T100_Tp40_SN80.dat" w l 
replot "bit_err_N160_K100_T100_Tp40_SN100.dat" w l 
