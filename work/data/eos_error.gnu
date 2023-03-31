reset

set view map

set form cb '10^{%L}'
set log cb

set xr [-8:0]
set yr [3:6]

set key right top opaque

xaxis='Q'
yaxis='T_from_eos_e'

set xl 'log(T/K)'
set yl 'Q'
pltfile='eostest_gasradrec.dat'
set multiplot layout 2,2

pl pltfile u (log10(column(xaxis))):(log10(column(yaxis))):\
   (abs(column('rel_error_for_eos_p_cf'))) w imag t 'Error on p for eos_p_cf'

pl pltfile u (log10(column(xaxis))):(log10(column(yaxis))):\
   (abs(column('rel_error_for_eos_p'))) w imag t 'Error on p for eos_p'

pl pltfile u (log10(column(xaxis))):(log10(column(yaxis))):\
   (abs(column('entropy_from_dp'))) w imag t 'Error on mu for eos_p'

#pl pltfile u (log10(column(xaxis))):(log10(column(yaxis))):\
   (abs(column('d_from_inventropy')/column('d')-1)) w imag t 'Error on mu for eos_p'

#pl pltfile u (log10(column(xaxis))):(log10(column(yaxis))):\
   (abs(column('mu_from_eos_p')-column('mu_from_eos_e'))) w imag t 'Error on mu for eos_p'

pl pltfile u (log10(column(xaxis))):(log10(column(yaxis))):\
   (abs(column('T_from_eos_p')-column('T_from_eos_e'))) w imag t 'Error on T for eos_p'

unset multiplot