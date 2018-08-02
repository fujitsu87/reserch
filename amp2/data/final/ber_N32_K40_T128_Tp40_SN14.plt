set logscale y
set xlabel "n"
set ylabel "BER"
plot "bit_err_N32_K40_T128_Tp40_SN14.dat" w lp linewidth 1 title "whole"
replot "bit_err_my_N32_K40_T128_Tp40_SN14.dat" w lp linewidth 1 title "my"
replot "bit_err_other_N32_K40_T128_Tp40_SN14.dat" w lp linewidth 1 title "other"
