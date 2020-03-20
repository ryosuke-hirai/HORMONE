reset
clear

set log y2
set y2tic
set ytic nomi

set yl 'v0 (km/s)'
set y2l 'delta M (Msun)'
set xl 'polar angle (degrees)'

#set term pngcairo
#set out 'parameter_distribution.png'

pl 'homology.dat' u 1:($3/1e5) w l lw 2 t 'v0',\
   '' u 1:2 w l axes x1y2 lw 2 t 'delta M'

