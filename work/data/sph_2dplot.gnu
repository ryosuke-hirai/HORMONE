reset

rsun = 6.96e10
msun = 1.989e33

set size ratio 1

box = 3000 #3e13/rsun
set xr [-box:box]
set yr [-box:box]

#set log cb
#set form cb '10^{%L}'

set pm3d map

rad = 'x1'
tht = 'x2'
val = 'o16'

do for [i=0:0:100]{

pltfile=sprintf('<paste gridfile.dat plt%8.8d.dat',i)
timelabel=sprintf('%d sec',i)
set label 1 timelabel at graph 0.5,graph 1.1 center
set style fill empty border rgb 'black'


    
splot \
pltfile \
   u ( column(rad)/rsun*sin(column(tht))):(column(rad)/rsun*cos(column(tht))):(column(val)) not w pm3d,\
'' u (-column(rad)/rsun*sin(column(tht))):(column(rad)/rsun*cos(column(tht))):(column(val)) not w pm3d

}