set logscale y
set xlabel "n"
set ylabel "BER"
plot "bit_err_N256_K32_T1000_Tp300_SN10.dat" w lp linewidth 1 title "whole"
replot "bit_err_my_N256_K32_T1000_Tp300_SN10.dat" w lp linewidth 1 title "my"
replot "bit_err_other_N256_K32_T1000_Tp300_SN10.dat" w lp linewidth 1 title "other"