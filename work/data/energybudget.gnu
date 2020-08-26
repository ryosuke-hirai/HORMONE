reset
clear
#set grid
set yr [-2:2]
set y2r [0:1e55]

set xl 'Time (ks)'
set yl 'Energy (10^{51}erg)'
set y2l 'Angular momentum (cgs)'

set y2tic
set ytic nomi
#set log y2

set term pngcairo size 12in,8in
set out 'energybudget.png'

pl 'angmom.dat'\
             u ($1/1e3):($9/1e51) w l lw 2 t 'Total energy',\
''       u ($1/1e3):($9*$12/1e51) w l lw 2 t 'Rotational energy',\
''          u ($1/1e3):($15/1e51) w l lw 2 t 'Internal energy',\
''          u ($1/1e3):($18/1e51) w l lw 2 t 'Gravitational energy',\
'' u ($1/1e3):(($21-$9*$12)/1e51) w l lw 2 t 'Kinetic energy (w/o rot)',\
'' u ($1/1e3):6 w l lw 2 dt '-' axes x1y2 t 'Angular momentum'