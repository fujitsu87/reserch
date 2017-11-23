set logscale y
set xlabel "n"
set ylabel "BER"
plot "bit_err_N256_K8_T1000_Tp64_SN10.dat" w l linewidth 4 title "whole"
replot "bit_err_my_N256_K8_T1000_Tp64_SN10.dat" w l linewidth 4 title "my"
replot "bit_err_other_N256_K8_T1000_Tp64_SN10.dat" w l linewidth 4 title "other"
