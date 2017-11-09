set logscale y
set xlabel "n"
set ylabel "BER"
plot "bit_err_N128_K64_T1024_Tp256_SN10.dat" w l linewidth 4 title "whole"
replot "bit_err_my_N128_K64_T1024_Tp256_SN10.dat" w l linewidth 4 title "my user"
replot "bit_err_other_N128_K64_T1024_Tp256_SN10.dat" w l linewidth 4 title "other user"