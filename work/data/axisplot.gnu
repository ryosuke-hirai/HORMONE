reset
clear

thres=1100000
ks=1000
set xr [0:2.2e8]
set cbr [0:(3e6-thres)/ks]
set palette defined (1 'yellow',2 'red', 3 'blue')

set term pngcairo size 15in,10in
set out 'similarity_axis.png'

unset colorbox
set multiplot layout 1,2
set xl 'r/t (cm/s)'
set yl 'Velocity (cm/s)'
pl '1Daveraged.data' ev :10 u ($2/($1-thres)):($1<3e6&&$1>thres?$4:NaN):(($1-thres)/ks) w l lw 3 lc palette not

set colorbox
set cbl 'Time (ks)'
set log y
set yr [*:1e14]
set yl 'Density*time^3 (g/cc*s^3)'
pl '1Daveraged.data' ev :10 u ($2/($1-thres)):($1<3e6&&$1>thres?($3*($1-thres)**3):NaN):(($1-thres)/ks) w l lw 3 lc palette not

unset multiplot