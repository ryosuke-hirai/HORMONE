reset

rsun = 6.96e10
msun = 1.989e33

set size ratio 1

box = 0.5 #3e13/rsun
set xr [-box:box]
set yr [-box:box]

#set log cb
#set form cb '10^{%L}'

set pm3d map

x1 = 'x1'
x2 = 'x2'
val = 'd'

do for [i=0:2200:100]{

pltfile(i)=sprintf('<paste gridfile.dat plt%11.11dms.dat',i)
timelabel=sprintf('%d sec',i)
set label 1 timelabel at graph 0.5,graph 1.1 center
set style fill empty border rgb 'black'

den='d'
splot \
pltfile(i) u (column(x1)):(column(x2)):(column(val)) w pm3d not

}