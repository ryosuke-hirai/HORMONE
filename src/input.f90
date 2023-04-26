module input_mod
 implicit none

 contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE READ_MESA
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To read MESA file

subroutine read_mesa(mesafile,r,m,rho,pres,comp,comp_list)

 use constants,only:msun,rsun

 implicit none

 character*100,intent(in):: mesafile
 character*30,allocatable::header(:),dum(:)
 real*8,allocatable,dimension(:),intent(out):: r, m, rho, pres
 real*8,allocatable,dimension(:,:),intent(out),optional:: comp
 character*10,allocatable,intent(out),optional::comp_list(:)
 real*8,allocatable,dimension(:,:):: dat
 character*10000 dumc
 integer nn,ui, i,j, lines, rows, nrel, nel,istat
 character(len=10):: element
 character(len=10),allocatable:: element_list(:)

 
!-----------------------------------------------------------------------------

! list of relevant elements ! ------------------------------------------------
 nrel = 20 ! number of elements in the list
 allocate(element_list(nrel))
 element_list = (/'h1','he3','he4','c12','n14','o16','ne20','mg24',&
                  'si28','s32','ar36','ca40','ti44','cr48','cr60',&
                  'fe52','fe54','fe56','co56','ni56'/)
 
! reading data from datafile ! -----------------------------------------------
 open(newunit=ui,file=mesafile,status='old',iostat=istat)
 if(istat/=0)then
  print*,'Error: Input MESA file cannot be found.'
  print*,'Make sure to specify the correct file name.'
  print'(2a)',"mesafile='",mesafile
  stop
 end if
 read(ui,'()')
 read(ui,'()')
 read(ui,*) lines, lines
 read(ui,'()')
 read(ui,'()')
 read(ui,'(a)') dumc

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

! find relevant chemical elements in the headers
 if(present(comp))then
  nel = 0
  do j = 1, nrel
   do i = 1, rows
    if(trim(header(i))==trim(element_list(j)))then
     nel = nel + 1 ! first count how many elements
     exit
    end if
   end do
  end do
  allocate(comp_list(nel))
  nn = 1
  element_loop: do j = 1, nrel
   do i = 1, rows
    if(trim(header(i))==trim(element_list(j)))then
     comp_list(nn) = header(i)
     nn = nn + 1
     if(nn>nel)exit element_loop
     exit
    end if
   end do
  end do element_loop
 end if

 do i = 1, lines
  read(ui,*) dat(lines-i+1,1:rows) 
 end do

 allocate(m(0:lines),r(0:lines),rho(0:lines),pres(0:lines),comp(1:nel,0:lines))
 do i = 1, rows
  if(trim(header(i))=='mass') m(1:lines) = dat(1:lines,i) * msun
  if(trim(header(i))=='density') rho(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='rho') rho(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='radius_cm') r(1:lines) = dat(1:lines,i)
  if(trim(header(i))=='radius') r(1:lines) = dat(1:lines,i) * rsun
  if(trim(header(i))=='logR') r(1:lines) = 10**dat(1:lines,i) * rsun
  if(trim(header(i))=='pressure') pres(1:lines) = dat(1:lines,i)
  if(present(comp))then
   do j = 1, nel
    if(trim(header(i))==trim(comp_list(j))) comp(j,1:lines) = dat(1:lines,i)
   end do
  end if
 end do

 close(ui)
 
 r(0) = 0d0
 m(0) = 0d0
 rho(0) = rho(1)
 pres(0) = pres(1)
 comp(:,0) = comp(:,1)

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

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                         SUBROUTINE ERROR_EXTRAS
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To output error message and stop simulation when failing to find extras file.

subroutine error_extras(simutype,extrasfile)

 character(len=*),intent(in):: simutype,extrasfile
 
!-----------------------------------------------------------------------------

 print*,'Error: Model parameter file cannot be found.'
 print'(5a)','Copy over "../para/extras_',simutype,'" to "',trim(extrasfile),&
             '" and specify model parameters'
 stop

end subroutine error_extras

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                          SUBROUTINE ERROR_NML
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To output error message and stop simulation when relevant namelist cannot be found in extras file

subroutine error_nml(simutype,extrasfile)

 character(len=*),intent(in)::simutype,extrasfile

!-----------------------------------------------------------------------------

 print*,'Error: extras file does not contain relevant namelist'
 print'(5a)','Copy over contents of "../para/extras_',simutype,'" to "',&
             trim(extrasfile),'" and specify model parameters'
 stop

return
end subroutine error_nml

end module input_mod
