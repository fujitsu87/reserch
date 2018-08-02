set logscale y
set xlabel "n"
set ylabel "BER"
plot "bit_err_N4_K2_T4_Tp2_SN1000000.dat" w lp linewidth 1 title "whole"
replot "bit_err_my_N4_K2_T4_Tp2_SN1000000.dat" w lp linewidth 1 title "my"
replot "bit_err_other_N4_K2_T4_Tp2_SN1000000.dat" w lp linewidth 1 title "other"
