reset
clear

istart=900000
iend=2220000
kizami=10000
rsun = 6.963e10
box=1500
set xr [-box:box]
set yr [-box:box]
set cbr [1e-14:1e-2]
set size square
set pm3d map
set log cb
set form cb '10^{%L}'
set palette define (0 'black',1 'blue',2 'green',3 'yellow',4 'orange',5 'red')

set term pngcairo size 15in,12in

do for [time=2030000:iend:kizami]{

    n = (time-istart)/1000
    pltfile=sprintf('<paste gridfile.dat plt%11.11ds.dat',time)
    outfile=sprintf('den%11.11ds.png',time)
    set label 1 sprintf('%dsec',time) at graph 0.5,graph 1.1 center
    set out outfile
    
    spl pltfile \
       u ($3/rsun) :($4/rsun):(abs($6)) not w pm3d,\
       '' u ($3/rsun):(-$4/rsun):($6*$12+$7!=0?(abs($6)):NaN) not w pm3d,\
       '' u (-$3/rsun):(-$4/rsun):($6*$12+$7!=0?(abs($6)):NaN) not w pm3d,\
       '' u (-$3/rsun):($4/rsun):(abs($6)) not w pm3d,\
    for [i=1:50]'particle_trajectories20Rsun.dat' ev ::n::n u (column(i*3)/rsun):(column(i*3+1)/rsun):(10) w p not lc 1 pt 7 ps 0.2,\
    for [i=1:50]'particle_trajectories30Rsun.dat' ev ::n::n u (column(i*3)/rsun):(column(i*3+1)/rsun):(10) w p not lc 2 pt 7 ps 0.2,\
    for [i=1:50]'particle_trajectories40Rsun.dat' ev ::n::n u (column(i*3)/rsun):(column(i*3+1)/rsun):(10) w p not lc 3 pt 7 ps 0.2,\
    for [i=1:50]'particle_trajectories50Rsun.dat' ev ::n::n u (column(i*3)/rsun):(column(i*3+1)/rsun):(10) w p not lc 4 pt 7 ps 0.2,\
    for [i=1:50]'particle_trajectories60Rsun.dat' ev ::n::n u (column(i*3)/rsun):(column(i*3+1)/rsun):(10) w p not lc 5 pt 7 ps 0.2,\
    for [i=1:50]'particle_trajectories20Rsun.dat' ev ::n::n u (-column(i*3)/rsun):(column(i*3+1)/rsun):(10) w p not lc 1 pt 7 ps 0.2,\
    for [i=1:50]'particle_trajectories30Rsun.dat' ev ::n::n u (-column(i*3)/rsun):(column(i*3+1)/rsun):(10) w p not lc 2 pt 7 ps 0.2,\
    for [i=1:50]'particle_trajectories40Rsun.dat' ev ::n::n u (-column(i*3)/rsun):(column(i*3+1)/rsun):(10) w p not lc 3 pt 7 ps 0.2,\
    for [i=1:50]'particle_trajectories50Rsun.dat' ev ::n::n u (-column(i*3)/rsun):(column(i*3+1)/rsun):(10) w p not lc 4 pt 7 ps 0.2,\
    for [i=1:50]'particle_trajectories60Rsun.dat' ev ::n::n u (-column(i*3)/rsun):(column(i*3+1)/rsun):(10) w p not lc 5 pt 7 ps 0.2,\
    for [i=1:50]'particle_trajectories20Rsun.dat' ev ::n::n u (column(i*3)/rsun):(-column(i*3+1)/rsun):(10) w p not lc 1 pt 7 ps 0.2,\
    for [i=1:50]'particle_trajectories30Rsun.dat' ev ::n::n u (column(i*3)/rsun):(-column(i*3+1)/rsun):(10) w p not lc 2 pt 7 ps 0.2,\
    for [i=1:50]'particle_trajectories40Rsun.dat' ev ::n::n u (column(i*3)/rsun):(-column(i*3+1)/rsun):(10) w p not lc 3 pt 7 ps 0.2,\
    for [i=1:50]'particle_trajectories50Rsun.dat' ev ::n::n u (column(i*3)/rsun):(-column(i*3+1)/rsun):(10) w p not lc 4 pt 7 ps 0.2,\
    for [i=1:50]'particle_trajectories60Rsun.dat' ev ::n::n u (column(i*3)/rsun):(-column(i*3+1)/rsun):(10) w p not lc 5 pt 7 ps 0.2,\
    for [i=1:50]'particle_trajectories20Rsun.dat' ev ::n::n u (-column(i*3)/rsun):(-column(i*3+1)/rsun):(10) w p not lc 1 pt 7 ps 0.2,\
    for [i=1:50]'particle_trajectories30Rsun.dat' ev ::n::n u (-column(i*3)/rsun):(-column(i*3+1)/rsun):(10) w p not lc 2 pt 7 ps 0.2,\
    for [i=1:50]'particle_trajectories40Rsun.dat' ev ::n::n u (-column(i*3)/rsun):(-column(i*3+1)/rsun):(10) w p not lc 3 pt 7 ps 0.2,\
    for [i=1:50]'particle_trajectories50Rsun.dat' ev ::n::n u (-column(i*3)/rsun):(-column(i*3+1)/rsun):(10) w p not lc 4 pt 7 ps 0.2,\
    for [i=1:50]'particle_trajectories60Rsun.dat' ev ::n::n u (-column(i*3)/rsun):(-column(i*3+1)/rsun):(10) w p not lc 5 pt 7 ps 0.2 	

    pause 0.01

}