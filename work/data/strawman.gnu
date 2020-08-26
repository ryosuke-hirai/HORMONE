reset
clear

set view map
set palette defined (0 'blue',1 'white',2 'red')
set palette define (0 'black',1 'blue',2 'green',3 'yellow',4 'orange',5 'red')

set log cb
rsun = 6.955e10
set form cb '{%L}'

set xr [0:2.1e1]
set yr [0:3000]

#set tmargin 0
#set bmargin 0
#set lmargin 0
#set rmargin 0

set cbr [1e-16:1e-2]#for density
#set cbr [-5e15:5e15]#for bernoulli
#set cbr [-5e7:5e7]#for velocity
set xl 'Time(10^5s)'
set term pngcairo size 18in,12in font 'Arial,8'
set out 'angulardep.png'

set multiplot layout 3,6

do for [i=0:85:5]{

    pltfile = sprintf('1Daveraged_%2.2d-%2.2ddeg.data',i,i+5)
    degrees = sprintf('%d-%ddegrees',i,i+5)
    
    set label 1 degrees at graph 0.5,graph 0.8 center front tc rgb 'white'
    
    
spl pltfile u (($1-9e5)/1e5):($3/rsun):4 w pm3d not

}

unset multiplot