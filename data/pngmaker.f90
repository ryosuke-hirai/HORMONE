implicit none

integer i,n
character*20 pngfile,datafile

n=500
   write(*,'(a)')"gnuplot << EOF"
!   write(*,'(a)')"set terminal png size 600,750"
   write(*,'(a)')"set terminal pngcairo size 15in, 11in enhanced"
!   write(*,'(a)')"set title 'logdensity plot'"
!   write(*,'(a)')"set view map"
!   write(*,'(a)')"set pm3d map"
   write(*,'(a)')"set logscale cb"
!   write(*,'(a)')"set logscale y"
!   write(*,'(a)')"set y2tics"
!   write(*,'(a)')"set ytics nomirror"
   write(*,'(a)')"set size ratio 0.666666666667"
!   write(*,'(a)')"set cbrange[1e-8:1e2]"!for density
   write(*,'(a)')"set cbrange[1e-7:1e0]"!for pressure
   write(*,'(a)')"set xrange[-12.6e11:12.6e11]"!for density
   write(*,'(a)')"set yrange[-8.4e11:8.4e11]"!for velocity
!   write(*,'(a)')"set palette rgbformulae 21,22,23"
   write(*,'(a)')"set palette define (0 'black',1 'blue',2 'green',3 'yellow',4 'orange',5 'red')"
!   write(*,'(a)')"set logscale y"
!   write(*,'(a)')"set format cb '%L'"
!   write(*,'(a)')"set noxtics"
!   write(*,'(a)')"set noytics"   
   write(*,'(a)')"set form cb '10^{%L}'"
!   write(*,'(a)')"set sample 1000000"   
!   write(*,'(a)')"box = 1.3e13"   
!   write(*,'(a)')"set xrange[-box:box]"
!   write(*,'(a)')"set xrange[-2e12:2e12]"
!   write(*,'(a)')"set yrange[-box/1.7:box/1.7]"
!   write(*,'(a)')"set yrange[1e12:1e20]"
!   write(*,'(a)')"mabiki = 1"
!   write(*,'(a)')'size = 0.3'
!   write(*,'(a)')'iro = 1'
!   write(*,'(a)')'kata = 7'
!   write(*,'(a)')'unset colorbox'
!   write(*,'(a)')"set xlabel ''"
do i=451,n
   write(pngfile,"(a,i7.7,a)")"den",i*100,".png"
   write(datafile,"(a,i7.7,a)")"plt",i*100,"s.dat"
   write(*,'(a,a,a)')"set output '",trim(pngfile),"'"
   write(*,'(a,a,a)')"plot '",trim(datafile),"' u 3:4:(\$5) t '' w image,''  u (-\$3):(\$4):(\$5) w image not"!,'' u (-\$3):4:(\$5) t '' w image,'' u (-\$3):(-\$4):(\$5) t '' w image" !for density
!   write(*,'(a,a,a)')"splot '",trim(datafile),"' ev :::2 u 3:4:(abs((\$11-\$10)/\$10)) t '' w image,'' ev:::2 u (\$3):(-\$4):(abs((\$11-\$10)/\$10)) w image not,'' ev :::2 u (-\$3):4:(abs((\$11-\$10)/\$10)) t '' w image,'' ev :::2 u (-\$3):(-\$4):(abs((\$11-\$10)/\$10)) t '' w image" !for density
!   write(*,'(5a)')"splot '< paste ",trim(datafile)," ../../headoncollision_final2/data/",trim(datafile),"' ev :::2 u 3:4:(abs(-\$10+\$20)/abs(\$10)) t '' w image,'' ev:::2 u (\$3):(-\$4):(abs(-\$10+\$20)/abs(\$10)) w image not,'' ev :::2 u (-\$3):4:(abs(-\$10+\$20)/abs(\$10)) t '' w image,'' ev :::2 u (-\$3):(-\$4):(abs(-\$10+\$20)/abs(\$10)) t '' w image" !for difference
!   write(*,'(a,a,a,a,a)')"plot '",trim(datafile),"' ev :::::1 u (-\$3):(\$7/\$5) lc 1 w l title'','' ev :::180 u 3:(\$7/\$5) lc 1 w l title '','plt0000000.dat' ev :::::1 u (-\$3):(\$7/\$5) lc 7 w l t '','' ev :::180 u 3:(\$7/\$5) w l lc 7 t ''" !for shock
!   write(*,'(a,a,a,a,a)')"pl '",trim(datafile),"' u (-\$3*cos(\$4)):(\$3*sin(\$4)):(((\$9)*sin(\$4)-(\$8)*cos(\$4))*1.5e4):(((\$8)*sin(\$4)+(\$9)*cos(\$4))*1.5e4) every 10:5 with vector lt 3 title '','",trim(datafile),"' u (-\$3*cos(\$4)):(-\$3*sin(\$4)):(((\$9)*sin(\$4)-(\$8)*cos(\$4))*1.5e4):(-((\$8)*sin(\$4)+(\$9)*cos(\$4))*1.5e4) every 10:5  with vector lt 3 title '', sqrt(4e13**2-x**2) notitle lt 1,-sqrt(4e13**2-x**2) notitle lt 1" !for velocity

!for particles
!!$   write(*,'(a)')'set multiplot'
!!$   write(*,'(a)')'set lmargin 0'
!!$   write(*,'(a)')'set rmargin 0'
!!$   write(*,'(a)')'set tmargin 0'
!!$   write(*,'(a)')'set bmargin 0'
!!$
!!$   write(*,'(a)')'set size noratio'
!!$   write(*,'(a)')'set xr [*:*]'
!!$   write(*,'(a)')'set yr [*:*]'   
!!$   write(*,'(a)')'set size 0.8,0.12'
!!$   write(*,'(a)')'set origin 0.1,0.8'
!!$   write(*,'(a)')'set ytics 1'
!!$   write(*,'(a)')'set xtics #offset 0,graph 0.1'
!!$   write(*,'(a)')'set xl "time (days)" #offset 0,1'
!!$
!!$   write(*,'(a,i4,a)')'pl "remfile" u ((\$2)/86400):3 w l t "remaining mass", "" u ((\$2)/86400):(\$1==',i,'000 ? \$3:1/0) w p pt 7 ps 1.5 t ""'
!!$
!!$
!!$   write(*,'(a)')'set size 0.8,1'
!!$   write(*,'(a)')'set origin 0.1,-0.1'
!!$   write(*,'(a)')'set size square'
!!$   write(*,'(a)')'set xr [-8e13:8e13]'
!!$   write(*,'(a)')'set yr [-8e13:8e13]'
!!$   write(*,'(a)')'set notics'
!!$   write(*,'(a)')'set xl ""'
!!$
!!$
!!$   write(*,'(a,a,a)')'plot "',trim(datafile),'" u (-(\$2)):(\$5>0 ? \$3:1/0) every mabiki pt kata ps size lt iro t "ejecta","" u (-(\$2)):(\$5<1 ? \$3:1/0) every mabiki pt kata ps size lt iro+1 t "star","" u (-(\$2)):(\$6>0 ? -(\$3):1/0) every mabiki pt kata ps size lt iro+2 t "bound","" u (-(\$2)):(\$6<1 ? -(\$3):1/0) every mabiki pt kata ps size lt iro+3 t "unbound"'
!!$
!!$   write(*,'(a)')'unset multiplot'
end do
   write(*,'(a)')"EOF"

end program
