set logscale y
set xlabel "Total number of iterations"
set ylabel "MSE"
plot "mse_h_N32_K40_T128_Tp40_SN8.dat" w l lc rgb "#0000ff" linewidth 3 notitle
replot "mse_h_N32_K40_T128_Tp40_SN10.dat" w l lc rgb "#0000ff" linewidth 3 notitle
replot "mse_h_N32_K40_T128_Tp40_SN11.dat" w l lc rgb "#0000ff" linewidth 3 notitle
replot "mse_h_N32_K40_T128_Tp40_SN12.dat" w l lc rgb "#0000ff" linewidth 3 notitle
replot "mse_h_N32_K40_T128_Tp40_SN14.dat" w l lc rgb "#0000ff" linewidth 3 notitle
replot "mse_h_N32_K40_T128_Tp40_SN16.dat" w l lc rgb "#0000ff" linewidth 3 notitle
replot "mse_h_N32_K40_T128_Tp40_SN18.dat" w l lc rgb "#0000ff" linewidth 3 notitle
replot "mse_h_N32_K40_T128_Tp40_SN20.dat" w l lc rgb "#0000ff" linewidth 3 notitle
replot "bit_err_N32_K40_T128_Tp40_SN8.dat" w l lc rgb "#ff0000" linewidth 3 notitle
replot "bit_err_N32_K40_T128_Tp40_SN10.dat" w l lc rgb "#ff0000" linewidth 3 notitle
replot "bit_err_N32_K40_T128_Tp40_SN11.dat" w l lc rgb "#ff0000" linewidth 3 notitle
replot "bit_err_N32_K40_T128_Tp40_SN12.dat" w l lc rgb "#ff0000" linewidth 3 notitle
replot "bit_err_N32_K40_T128_Tp40_SN14.dat" w l lc rgb "#ff0000" linewidth 3 notitle
replot "bit_err_N32_K40_T128_Tp40_SN16.dat" w l lc rgb "#ff0000" linewidth 3 notitle
replot "bit_err_N32_K40_T128_Tp40_SN18.dat" w l lc rgb "#ff0000" linewidth 3 notitle
replot "bit_err_N32_K40_T128_Tp40_SN20.dat" w l lc rgb "#ff0000" linewidth 3 notitle