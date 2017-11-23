set logscale y
set xlabel "n"
set ylabel "BER"
plot "bit_err_N128_K32_T128_Tp32_SN10.dat" w lp linewidth 4 title "whole"
replot "bit_err_my_N128_K32_T128_Tp32_SN10.dat" w lp linewidth 4 title "my"
replot "bit_err_other_N128_K32_T128_Tp32_SN10.dat" w lp linewidth 4 title "other"
