reset
clear

set xr [10000:*]
set xr [0.1:*]
set yr [0:14]
#set log x

set term pngcairo
set out 'Unboundmasses.png'

set yl 'Unbound mass (Msun)'
set xl 'Time (s)'

pl 'angmom.dat' u ($1-9e5):(100-$3) w l lw 2 lc 1 t 'Tinj=0.2orb','' u ($1-9e5):(100-$4) w l dt '-' lw 2 lc 1 not,\
   '../../aeff40_x05_tinj0.4_Eheat1.0/data/angmom.dat' u ($1-9e5):(100-$3) w l lw 2 lc 2 t 'Tinj=0.4orb','' u ($1-9e5):(100-$4) w l dt '-' lw 2 lc 2 not,\
   '../../aeff40_x05_tinj1_Eheat1.0/data/angmom.dat' u ($1-9e5):(100-$3) w l lw 2 lc 3 t 'Tinj=1orb','' u ($1-9e5):(100-$4) w l dt '-' lw 2 lc 3 not,\
   '../../aeff40_x05_tinj5_Eheat1.0/data/angmom.dat' u ($1-9e5):(100-$3) w l lw 2 lc 4 t 'Tinj=5orb','' u ($1-9e5):(100-$4) w l dt '-' lw 2 lc 4 not,15 w l lw 2 lc rgb 'black' t 'Total energy',15 w l dt '-' lw 2 lc rgb 'black' t 'Bernoulli'