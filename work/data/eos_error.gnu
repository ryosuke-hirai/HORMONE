reset

set view map

set form cb '10^{%L}'
set log cb

set xr [3:6]
set yr [-8:0]

set key right top opaque

xaxis='T_from_eos_e'
yaxis='Q'

set xl 'log(T/K)'
set yl 'Q'

set multiplot layout 2,2

pl 'eos.dat' u (log10(column(xaxis))):(log10(column(yaxis))):\
   (abs(column('rel_error_for_eos_p_cf'))) w imag t 'Error on p for eos_p_cf'

pl 'eos.dat' u (log10(column(xaxis))):(log10(column(yaxis))):\
   (abs(column('rel_error_for_eos_p'))) w imag t 'Error on p for eos_p'

pl 'eos.dat' u (log10(column(xaxis))):(log10(column(yaxis))):\
   (abs(column('mu_from_eos_p')-column('mu_from_eos_e'))) w imag t 'Error on mu for eos_p'

pl 'eos.dat' u (log10(column(xaxis))):(log10(column(yaxis))):\
   (abs(column('T_from_eos_p')-column('T_from_eos_e'))) w imag t 'Error on T for eos_p'

unset multiplot