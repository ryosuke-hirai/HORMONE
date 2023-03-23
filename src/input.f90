module input_mod
 implicit none

 contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE READ_MESA
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To read MESA file

subroutine read_mesa(mesafile,r,m,rho,pres)

 use constants,only:msun

 implicit none

 character*100,intent(in):: mesafile
 character*30,allocatable::header(:),dum(:)
 real*8,allocatable,dimension(:),intent(out):: r, m, rho, pres
 character*10000 dumc
 integer nn, nol, ii, i, rows
 
!-----------------------------------------------------------------------------

 ! reading data from datafile ! -----------------------------------------------
 open(unit=40,file=mesafile,status='old')
 read(40,'()')
 read(40,'()')
 read(40,*) lines, lines
 read(40,'()')
 read(40,'()')
 read(40,'(a)') dumc

! counting rows
 allocate(dum(500)) ; dum = 'aaa'
 read(dumc,*,end=101) dum
101 do i = 1, 500
  if(dum(i)=='aaa')then
   rows = i - 1
   exit
  end if
 end do

 allocate(header(1:rows),dat(1:lines,1:rows))
 header(1:rows) = dum(1:rows)
 deallocate(dum)

 do i = 1, lines
  read(40,*) dat(lines-i+1,1:rows) 
 end do

 allocate(m(0:lines),r(0:lines),rho(0:lines),pres(1:lines))
 do i = 1, rows
  if(trim(header(i))=='mass') m(1:lines) = dat(1:lines,i) * msun
  if(trim(header(i))=='density') rho(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='radius_cm') r(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='pressure') pres(1:lines) = dat(1:lines,i)
 end do

 r(0) = 0d0
 m(0) = 0d0
 rho(0) = rho(1)

return
end subroutine read_mesa

! set ejecta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!!$ ejtbinfile = 'ejtbin17_5.0B_dis1.62E+15cm.dat'
!!$ sep = 4000d0*rsun
!!$
!!$ open(unit=90, file=ejtbinfile, status='old',form='unformatted')
!!$ count = -2
!!$ do
!!$  read(90,end=311)
!!$  count = count+1
!!$ end do
!!$311 close(90)!counting lines of ejtfile
!!$
!!$ ejdatnum = count
!!$ allocate( t_ej(count),d_ej(count),p_ej(count),e_ej(count), &
!!$           v_ej(count),m_ej(count) )
!!$
!!$ open(unit=90, file=ejtbinfile, status='old',form='unformatted')
!!$
!!$ read(90)psmass
!!$ read(90)ejectadistance
!!$
!!$ do i = 1, count
!!$  read(90)t_ej(i), d_ej(i), p_ej(i), e_ej(i), v_ej(i), m_ej(i)
!!$  if(tstart==0.d0.and.v_ej(i)>=1.d4)then
!!$  !if(tstart==0d0.and.d_ej(i)>6d-16)then
!!$   tstart  = t_ej(i)
!!$   tstartn = i
!!$  end if
!!$  if(p_ej(i)>pmax) pmax = p_ej(i)
!!$ end do
!!$ close(90) !reading ejecta data
!!$
!!$ allocate( nsdis(is-2:ie+2,js-2:je+2,ks-2:ke+2), &
!!$           nsdfr(is-2:ie+2,js-2:je+2,ks-2:ke+2), &
!!$           nssin(is-2:ie+2,js-2:je+2,ks-2:ke+2), &
!!$           nscos(is-2:ie+2,js-2:je+2,ks-2:ke+2) )
!!$
!!$ j = js
!!$ do k = ks-2,ke+2
!!$  do i = is-2,ie+2
!!$   nsdis(i,j,k) = sqrt(x1(i)*x1(i)+(sep-x3(k))**2d0)
!!$   nsdfr(i,j,k) = ejectadistance / nsdis(i,j,k)
!!$   nscos(i,j,k) = (sep-x3(k)) / nsdis(i,j,k)
!!$   nssin(i,j,k) = x1(i) / nsdis(i,j,k)
!!$  end do
!!$ end do


end module input_mod
