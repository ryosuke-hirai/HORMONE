reset
clear

set style fill solid

#set term pngcairo
#set out 'mass_distribution_histogram.png'

#pl \
   'latitude_distribution_v100.dat' u ($1-2.5):2:(5) w boxes lc pal t 'v>100km/s',\
   'latitude_distribution_v200.dat' u ($1-2.5):2:(4) w boxes lc pal t 'v>200km/s',\
   'latitude_distribution_v300.dat' u ($1-2.5):2:(3) w boxes lc pal t 'v>300km/s',\
   'latitude_distribution_v400.dat' u ($1-2.5):2:(2) w boxes lc pal t 'v>400km/s',\
   'latitude_distribution_v500.dat' u ($1-2.5):2:(1) w boxes lc pal t 'v>500km/s'

unset colorbox

set terminal gif size 640,640 animate delay 3 
set out 'angle_timeevo.gif'
set lmargin 0
set rmargin 0
set bmargin 0
set tmargin 0

inifile = 1000000

do for [n=inifile:inifile+1660000:10000]{

 set multiplot

set xr [0:1]
set yr [0:1.5]
set y2r [0:1.4e48]
#set cbr [0:5]
set xl 'cos(theta)'
set yl 'Unbound mass (Msun)'
set y2l 'Kinetic energy (erg)'
set y2tic
set ytic nomi
set key left top
    
 set origin 0.2,0.1
 set size 0.6,0.55
    
    m=(n-inifile)/1000

 set label 1 sprintf('%dsec',n-inifile) at graph 0.9, graph 0.95 right
 set arrow 1 from sqrt(2)/2, graph 0 to sqrt(2)/2, graph 1 nohead dt '-' front
 set arrow 2 from sqrt(3)/2, graph 0 to sqrt(3)/2, graph 1 nohead dt '-' front

 set label 2 '45°' at sqrt(2)/2,graph 0.8 rotate by 90
 set label 3 '30°' at sqrt(3)/2,graph 0.8 rotate by 90
    
pl 'costheta_distribution.dat' ev :::m::m u 1:2 w boxes t 'v>0km/s' lc 1,\
   '' ev :::m::m u 1:3 w boxes t 'v>100km/s' lc 2,\
   '' ev :::m::m u 1:4 w boxes t 'v>200km/s' lc 3,\
   '' ev :::m::m u 1:5 w boxes t 'v>300km/s' lc 4,\
   '' ev :::m::m u 1:6 w boxes t 'v>400km/s' lc 5,\
   '' ev :::m::m u 1:7 w boxes t 'v>500km/s' lc 6,\
   '' ev :::m::m u 1:8 w l axes x1y2 t 'Kinetic energy' lc rgb 'black' lw 2

 set origin 0.2,0.75
 set size 0.6,0.18

 unset label 1
 unset arrow 1
 unset arrow 2
 unset label 2
 unset label 3

set xr [0:1700]
set yr [0:14]
unset y2tic
set xl 'Time (ks)'
set yl 'Unbound mass (Msun)'
set y2l ''

pl 'angmom.dat' u ($1/1e3):(100-$3) w l not,'' u ($1==n-inifile?$1/1e3:NaN):(100-$3) w p ps 2 pt 7 not

 unset multiplot
    
    pause 0.02
}


set out