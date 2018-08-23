set logscale y
set xlabel "n"
set ylabel "BER"
plot "bit_err_N8_K32_T128_Tp32_SN10.dat" w lp linewidth 1 title "whole"
replot "bit_err_my_N8_K32_T128_Tp32_SN10.dat" w lp linewidth 1 title "my"
replot "bit_err_other_N8_K32_T128_Tp32_SN10.dat" w lp linewidth 1 title "other"
