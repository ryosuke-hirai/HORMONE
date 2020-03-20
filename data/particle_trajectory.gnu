reset
clear

box=1000
set xr [0:box]
set yr [0:box]
set size square
set palette rgb 33,13,10
rsun= 6.963e10

set xl 'r (Rsun)'
set yl 'z (Rsun)'
set cbl 'time (sec)'

#set term pngcairo

#set out 'trajectory_20Rsun.png'
#set label 1 'r=20Rsun' at graph 0.95, graph 0.9 right
#pl for [i=1:50:1]'particle_trajectories20Rsun.dat' u (column(i*3)/rsun):(column(i*3+1)/rsun):1 w l not lc pal
#set out 'trajectory_30Rsun_zoomin.png'
#set label 1 'r=30Rsun' at graph 0.95, graph 0.9 right
#pl for [i=1:50:1]'particle_trajectories30Rsun.dat' u (column(i*3)/rsun):(column(i*3+1)/rsun):1 w l not lc pal
#set out 'trajectory_40Rsun.png'
set label 1 'r=40Rsun' at graph 0.95, graph 0.9 right
pl for [i=1:50:1]'particle_trajectories40Rsun.dat' u (column(i*3)/rsun):(column(i*3+1)/rsun):1 w l not lc pal
#set out 'trajectory_50Rsun.png'
#set label 1 'r=50Rsun' at graph 0.95, graph 0.9 right
#pl for [i=1:50:1]'particle_trajectories50Rsun.dat' u (column(i*3)/rsun):(column(i*3+1)/rsun):1 w l not lc pal
#set out 'trajectory_60Rsun.png'
#set label 1 'r=60Rsun' at graph 0.95, graph 0.9 right
#pl for [i=1:50:1]'particle_trajectories60Rsun.dat' u (column(i*3)/rsun):(column(i*3+1)/rsun):1 w l not lc pal
#set out 'trajectory_70Rsun.png'
#set label 1 'r=70Rsun' at graph 0.95, graph 0.9 right
#pl for [i=1:50:1]'particle_trajectories70Rsun.dat' u (column(i*3)/rsun):(column(i*3+1)/rsun):1 w l not lc pal