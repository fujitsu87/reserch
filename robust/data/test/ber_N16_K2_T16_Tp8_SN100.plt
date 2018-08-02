set logscale y
set xlabel "n"
set ylabel "BER"
plot "bit_err_N16_K2_T16_Tp8_SN100.dat" w lp linewidth 1 title "whole"
replot "bit_err_my_N16_K2_T16_Tp8_SN100.dat" w lp linewidth 1 title "my"
replot "bit_err_other_N16_K2_T16_Tp8_SN100.dat" w lp linewidth 1 title "other"
