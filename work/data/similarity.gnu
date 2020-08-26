reset
clear

set log y

initime=7e5
lastime=1.579e6
kizami=4e4

set xr [0:4e8]
set yr [1:1e14]
set form y '10^{%L}'
set xl 'r/t (cm/s)'
set yl 'rho*t^3 (cgs)'
#f(x)=a*exp(-x/b)
f(x)=-x/b+lna
b=1e7
a=2e34/8./pi/b**3
lna=log(a)
msun=1.986e33

#set term pngcairo

#set multiplot layout 5,2
set print 'homology.dat'
do for [j=1:200]{
#fit [0:1.e8][1:1e16] f(x) sprintf('<paste gridfile.dat plt%11.11ds.dat',lastime) ev :::j::j u ($3/(lastime-initime)):($6*$12+$7>0?$6*(lastime-initime)**3:NaN) via a,b
fit [0:1e8][1:1e16] f(x) sprintf('<paste gridfile2.dat plt%11.11ds.dat',lastime) ev :::j::j u ($3/(lastime-initime)):($6*$12+$7>0?log($6*(lastime-initime)**3):NaN) via lna,b
a=exp(lna)
    set label 1 sprintf('theta= %.2f deg',90/pi*(acos(1-j/200.)+acos(1-(j-1)/200.))) at graph 0.5,graph 0.9 center
    set label 2 sprintf('\delta M=%.2f Msun',a*8*pi*b**3/msun) at graph 0.5,graph 0.85 center
    set label 3 sprintf('v_0=%.2f km/s',b/1e5) at graph 0.5, graph 0.8 center
#    if (j>120){initime=1e6}
#    set out sprintf('ejecta_fit%3.3d.png',j)
print  90/pi*(acos(1-j/200.)+acos(1-(j-1)/200.)),a*8*pi*b**3/msun, b
#   pl for [i=initime+kizami*10:lastime:kizami] sprintf('<paste gridfile2.dat plt%11.11ds.dat',i) ev :::j::j u ($3/(i-initime)):($6*$12+$7>=0?$6*(i-initime)**3:NaN) w l lc i/kizami lw 2 not,\
    for [i=initime+kizami*10:lastime:kizami] sprintf('<paste gridfile2.dat plt%11.11ds.dat',i) ev :::j::j u ($3/(i-initime)):($6*(i-initime)**3) w l lc i/kizami t sprintf('%dks',(i-initime)/1e3),\
    a*exp(-x/b) lc rgb 'black' dt '-' lw 2 t 'fit'

}

#unset multiplot
