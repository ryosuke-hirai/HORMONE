reset
clear

set size ratio 1
rsun =6.96e10
box=700
set xr [-box:box]
set yr [-box:box]
set yr [0:2*box]
#set term pngcairo# gif animate delay 100
#set out 'ptcmovie.gif'

do for [n=2400000:2660000:10000]{

ptcfile=sprintf('ptc%11.11ds.dat',n)
outfile=sprintf('ptc%11.11ds.png',n)
timelabel=sprintf('%d sec',n)
set label 1 timelabel at graph 0.5,graph 0.95 center front
#set out outfile
#plot ptcfile \
   ev 5 u ( $5):($3==1?$6:1/0) w dots not,'' ev 5 u ( $5):($3==0?$6:1/0) w dots not,\
'' ev 5 u (-$5):($3==1?$6:1/0) w dots not,'' ev 5 u (-$5):($3==0?$6:1/0) w dots not
pn1=200
pn2=1100
pn3=35150
plot ptcfile \
   ev 1 u ( $5/rsun):( $6/rsun):(-$3+6) w dots not lc var,\
'' ev 1 u ( $5/rsun):(-$6/rsun):(-$3+6) w dots not lc var,\
'' ev 1 u (-$5/rsun):( $6/rsun):(-$3+6) w dots not lc var,\
'' ev 1 u (-$5/rsun):(-$6/rsun):(-$3+6) w dots not lc var,\
'' ev ::pn1::pn1 u ($5/rsun):($6/rsun) w p not pt 7 lc rgb 'black',\
'' ev ::pn3::pn3 u ($5/rsun):($6/rsun) w p not pt 7 lc rgb 'black',\
'' ev ::pn2::pn2 u ($5/rsun):($6/rsun) w p not pt 7 lc rgb 'red' ,x lc rgb 'grey' not,0.4*x

}
