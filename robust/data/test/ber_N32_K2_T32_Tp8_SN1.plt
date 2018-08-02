set logscale y
set xlabel "n"
set ylabel "BER"
plot "bit_err_N32_K2_T32_Tp8_SN1.dat" w lp linewidth 1 title "whole"
replot "bit_err_my_N32_K2_T32_Tp8_SN1.dat" w lp linewidth 1 title "my"
replot "bit_err_other_N32_K2_T32_Tp8_SN1.dat" w lp linewidth 1 title "other"
