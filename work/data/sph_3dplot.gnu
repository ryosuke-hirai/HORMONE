reset

rsun = 6.96e10
msun = 1.989e33

set size ratio 1

box = 70 #3e13/rsun
amp = 0e-7

set xr [-box:box]
set yr [-box:box]
set zr [0:box]

set log cb
#set form cb '10^{%L}'

rad = 'x1'
tht = 'x2'
phi = 'x3'
v1  = 'v1'
v3  = 'v3'
val = 'd'
#set palette defined (-1 'blue',0 'white',1 'red')
#load '/mnt/hgfs/shared/palette/inferno.pal'
#set cbr [2.3e13:2.35e13]
#unset log cb
set cbr [1e-10:2e0]
#set cbr [0:2e7]
#set term gif anim size 640,640
#set out 'movie.gif'
set form cb '10^{%L}'
set cbl 'Density (g/cc)'
#set term pngcairo size 640,640
#set pm3d map
set xyplane 0
set view equal xyz
set view 45,45

set style fill empty border rgb 'black'
do for [i=00:00:100]{

    pltfile=sprintf('<paste gridfile.dat plt%11.11dmin.dat',i)
    othfile=sprintf('<paste othergrid.dat oth%11.11dmin.dat',i)
    verfile=sprintf('<paste othergrid.dat ver%11.11dmin.dat',i)
    pngfile=sprintf('png/den3d%11.11dmin.png',i)
#set out pngfile
sinkfile='<tail -n +2 sinks.dat'
timelabel=sprintf('%d min',i)
set label 1 timelabel at screen 0.5,screen 0.95 center
set style fill empty border rgb 'black'

set key autotitle columnhead

sinkx(i,n)=sprintf('sink_%d_x%d',n,i)
    
#splot \
pltfile \
   u ( column(rad)/rsun*cos(column(phi))):(column(rad)/rsun*sin(column(phi))):(0):(abs(column(val))) not w pm3d,\
   for [n=1:2] sinkfile u (abs(column('time')-i*60)<1e-3?column(sinkx(1,n))/rsun:NaN):(column(sinkx(2,n))/rsun):(0) w p pt 7 ps 0.5 lc rgb 'blue'#,\
pltfile ev 30:30 u ( column(rad)/rsun*cos(column(phi))):(column(rad)/rsun*sin(column(phi))):(0):((column(v1)*cos(column(phi))-column(v3)*sin(column(phi)))*amp):((column(v1)*sin(column(phi))+column(v3)*cos(column(phi)))*amp):(0) w vec lc rgb 'white'#,\
othfile \
	  u (column(rad)/rsun*cos(column(tht))):(column(rad)/rsun*sin(column(tht))):(abs(column(val))) not w pm3d

# Vertical slice
#    splot \
	  othfile \
	  u (column(rad)/rsun*sin(column(tht))):(column(rad)/rsun*column(cos(tht))):(abs(column(val))) not w pm3d

    # Three-pane plot
    set arrow 1 from -box,0,0 to box,0,0 lc rgb 'white' nohead front
    set arrow 2 from 0,-box,0 to 0,box,0 lc rgb 'white' nohead front
    splot \
	  pltfile u (column(rad)/rsun*cos(column(phi))):(column(rad)/rsun*sin(column(phi))):(0):(column(val)) not w pm3d,\
	  othfile u (column(rad)/rsun*sin(column(tht))):(box):(column(rad)/rsun*cos(column(tht))):(column(val)) not w pm3d,\
	  verfile u (-box):(column(rad)/rsun*sin(column(tht))):(column(rad)/rsun*cos(column(tht))):(column(val)) not w pm3d,\
	  for [n=2:2] sinkfile u (abs(column('time')-i*60)<1e-3?column(sinkx(1,n))/rsun:NaN):(column(sinkx(2,n))/rsun):(0) w p pt 7 ps 0.5 lc rgb 'blue'#,\

print i
}
#set out
