reset
clear

box = 24000
hi = 1.66666666667
set xr [-box/hi:box/hi]
set yr [-box:box]
set size ratio hi

set xl 'r (au)'
set yl 'z (au)'

load 'viridis.pal'

#set term pngcairo transparent size 480, 640
#set out '1dhomunculus_1.5e-2wind.png'
#set cbr [0:25]
set log cb
set cbr [1:50]
pl 'analytic_homunculus.dat' \
   u ($1<90?( $2*sin($1/180*pi)):NaN):( $2*cos($1/180*pi)):4 w l lc pal lw 3 not,\
'' u ($1<90?(-$2*sin($1/180*pi)):NaN):( $2*cos($1/180*pi)):4 w l lc pal lw 3 not,\
'' u ($1<90?( $2*sin($1/180*pi)):NaN):(-$2*cos($1/180*pi)):4 w l lc pal lw 3 not,\
'' u ($1<90?(-$2*sin($1/180*pi)):NaN):(-$2*cos($1/180*pi)):4 w l lc pal lw 3 not,\
   'stan_homunculus.dat' \
   u ($1<90?( $2*sin($1/180*pi)):NaN):( $2*cos($1/180*pi)):4 w l lc pal lw 3 not,\
'' u ($1<90?(-$2*sin($1/180*pi)):NaN):( $2*cos($1/180*pi)):4 w l lc pal lw 3 not,\
'' u ($1<90?( $2*sin($1/180*pi)):NaN):(-$2*cos($1/180*pi)):4 w l lc pal lw 3 not,\
'' u ($1<90?(-$2*sin($1/180*pi)):NaN):(-$2*cos($1/180*pi)):4 w l lc pal lw 3 not

