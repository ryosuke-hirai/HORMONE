reset
#clear

box = 1
set xr [0:box]
set yr [0:box]
set cbr [0.01:0.4]
set size square

#set log cb

pltfile(i)=sprintf('<paste gridfile.dat plt%11.11dms.dat',i)

do for [i=0:500:10]{

timelabel=sprintf('%d ms',i)
set label 1 timelabel at graph 0.5,graph 0.95 center front

plot pltfile(i) u (column('x1')):(column('x3')):(column('d')) w imag

}
