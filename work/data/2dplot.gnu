reset
#clear

rsun = 6.96e10
msun = 1.989e33
inifile = 700000
set size ratio 4./4 #0.66666666667#0.5 #333333333
set palette rgbformulae 21,22,23
set palette define (0 'black',1 'blue',2 'green',3 'yellow',4 'orange',5 'red')
#set palette define (0 'blue', 1 'white', 2 'red')
#set palette maxcolors 10
box = 300 #3e13/rsun
set xr [-box:box]
set yr [-box:box]

set log cb
set form cb '10^{%L}'

#set term pngcairo size 15in,12in
#set term gif animate delay 50 optimize size 640,480
#set out 'test1_density.gif'
#set out 'tinj1_vel.gif'
#set out 'tinj01_closeup.gif'
set pm3d map

do for [i=809000:809000:10000]{

pltfile=sprintf('<paste gridfile2.dat plt%11.11ds.dat',i)
ptcfile=sprintf('ptc%11.11ds.dat',i)
outfile=sprintf('vel%11.11ds.png',i)
timelabel=sprintf('%d sec',i-inifile)
set label 1 timelabel at graph 0.5,graph 1.1 center
set style fill empty border rgb 'black'
    amp=1e5/rsun

set contour base 
year=24*3600*365.25
Mdot=1e-3*msun/year
vinf=6e7
amu=1.6605402e-24
arad=7.5646e-15
kbol=1.380658e-16
h = 6.6260755e-27
set cntrparam levels discrete 0
set linetype 2 lc rgb 'white'

set cntrparam levels discrete 0 #1e4,1e5,1e6,1e7,1e8,1e9,1e10
set cntrlabel onecolor

#    ((log((2.*pi*amu*kbol/h**2.)**1.5*$13**1.5/$6*$14)+2.5-5.95e1)/$14 +4.*arad/3.*$13**3./$6/kbol*amu)entropy
#set cbr [0:1]
#set cbr [15:60]
#    set cbr [0:5e-5]
#set object 1 circle at 0,0 size 500 front
#set cbr [100:10000]
splot \
pltfile \
   u ($3/rsun*sin($4)) :($3/rsun*cos($4)):(0):($6*$12+$7!=0.1?(($6)):NaN) not w pm3d,\
'' u (-$3/rsun*sin($4)):($3/rsun*cos($4)):(0):($6*$12+$7<0.?(($6)):NaN) not w pm3d,\
'' u ($3/rsun*sin($4)):(-$3/rsun*cos($4)):(0):($6*$12+$7!=0.1?(($6)):NaN) not w pm3d,\
'' u (-$3/rsun*sin($4)):(-$3/rsun*cos($4)):(0):(($6)) not w pm3d#,\
'' u ($3/rsun*sin($4)) :($3/rsun*cos($4)):(($9**2+$10**2+$11**2)/(2*$12)+1) not w pm3d nosurface lw 2 lc 2 #,\
    '' u ($3/rsun*sin($4)):($3/rsun*cos($4)):(0):(($9*sin($4)+$10*cos($4))*amp):(($9*cos($4)-$10*sin($4))*amp):(1) ev 10:2:100:1::7 not w vector lc rgb 'white',\
'' u ($3/rsun*sin($4)):($3/rsun*cos($4)):(0):(($9*sin($4)+$10*cos($4))*amp):(($9*cos($4)-$10*sin($4))*amp):(1) ev 10:10:100:7 not w vector lc rgb 'white'




Mtot=70*1.989e33
fac=0.5*pi*(7*6.96e10)**2.
G=6.7e-8
#set contour base
#set style textbox opaque margins  0.5,  0.5 fc  bgnd noborder# lw  1.0
#set cntrparam levels discrete 0.5, 1.0, 10.0, 100.0, 1000
#set cntrlabel onecolor
#set cntrlabel  #format '%8.3g' font ',7' start 5 interval 20
#set cntrlabel font ",7" front
#unset surface
#splot \
pltfile \
   u ($3/rsun) :($4/rsun):(abs(sqrt(G**3.*Mtot**5.)*0.25/(2*pi*($3*2)**2.5*$6*($10-sqrt(-0.5*$12))**3.*fac))) not w pm3d lw 2,\
'' u ($3/rsun):(-$4/rsun):($6*$12+$7!=0?(abs(sqrt(G**3.*Mtot**5.)*0.25/(2*pi*($3*2)**2.5*$6*($10-sqrt(-0.5*$12))**3.*fac))):NaN) not w pm3d lw 2,\
'' u (-$3/rsun):(-$4/rsun):($6*$12+$7!=0?(abs(sqrt(G**3.*Mtot**5.)*0.25/(2*pi*($3*2)**2.5*$6*($10-sqrt(-0.5*$12))**3.*fac))):NaN) not w pm3d lw 2,\
'' u (-$3/rsun):($4/rsun):(abs(sqrt(G**3.*Mtot**5.)*0.25/(2*pi*($3*2)**2.5*$6*($10-sqrt(-0.5*$12))**3.*fac))) not w pm3d lw 2
    
}

set out