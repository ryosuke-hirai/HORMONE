reset
#clear
set yr [*:*]
set log y
#set log y2
#set log x
#set y2tics

set grid

#set form y '{%L}'
#set xr [0:100]
#set xr [1.11e8:5e14]

#set log

rsun = 6.96e10
msun = 1.989e33

do for [t=000:1000000:10000] {
#  outfile = sprintf('animation/bessel%03.0f.png',t)
#  set output outfile
  pfile = sprintf('<paste gridfile.dat plt%11.11ds.dat',t)
  plot '<paste gridfile.dat plt00000000000s.dat' u ($2/rsun):($4)  w l,pfile u ($2/rsun):($4) w l lw 2#,'' u ($2/6.96e10):($7) w l axes x1y2#,\
'' u ($2/rsun):7 w l lw 2 axes x1y2,0 lw 2 axes x1y2
       
#plot '<paste gridfile.dat plt00000000000s.dat' u ($11/1.989e33):4 w l,pfile u ($11/1.989e33):4 w l lw 2,\
       '../100Msunprofile.data' ev ::7 u 2:18 w l #axes x1y2
       
  pause 0.1
  }

#pl '<paste gridfile.dat plt00000003000s.dat' u 2:7 w lp#,x**(-3)*4e33#,\
'../100Msunprofile.data' ev ::7 u ($12*6.96e10):20 w l
