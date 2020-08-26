reset
clear

r =5
omega=2
vr=1
end=900
kizami=1

set sample 10000
set size ratio 0.5
set xr [-10:10]
set yr [0:10]

do for [i=1:end:kizami]{

theta=i/1800.*pi
    
vphi=omega*sin(theta)
#vphi=1000*(r*sin(theta))**1.5

set arrow i from r*sin(theta),r*cos(theta) to r*sin(theta)+sqrt(vphi**2+vr*vr*sin(theta)**2),r*cos(theta)+vr*cos(theta)

}

do for [i=1:end:kizami]{

theta=i/1800.*pi
    
vphi=0

set arrow i+end from r*sin(theta),r*cos(theta) to r*sin(theta)+sqrt(vphi**2+vr*vr*sin(theta)**2),r*cos(theta)+vr*cos(theta) lc rgb 'red'

}

pl sqrt(1-x**2)