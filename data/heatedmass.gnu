reset
clear

set term pngcairo
set out 'heatedmass.png'

set xl 'Time (s)'
set yl 'M(r<15Rsun) (Msun)'

pl 'timeevo.dat' u ($1-9e5):($9*2) w l lc 1 not,\
   '' u ($1<=924000?($1-9e5):NaN):($9*2) w l lw 5 lc 1 t 'Tinj=0.2orb',\
   '../../aeff40_x05_tinj0.4_Eheat1.0/data/timeevo.dat' u ($1-9e5):($9*2) w l lc 2 not,\
   '' u ($1<=949000?($1-9e5):NaN):($9*2) w l lw 5 lc 2 t 'Tinj=0.4orb'