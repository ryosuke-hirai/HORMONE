!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE EXPLOSION
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set initial condition for 1D explosion.

subroutine explosion

 use grid
 use physval
 use constants
 use gravmod,only:mc

 implicit none

 integer lines, rows, n
 character*10000 dumc
 character*20,allocatable:: header(:),dum(:)
 real*8,allocatable,dimension(:,:):: dat
 real*8,allocatable,dimension(:):: m, r, pres
 real*8 Eexp, PNSmass, Ebind

!-----------------------------------------------------------------------------

! setting parameters for explosion ! -----------------------------------------
 Eexp = 1.d51 ! erg
 PNSmass = 1.4d0 * msun

! reading data from datafile ! -----------------------------------------------
 open(unit=40,file='progenitor025.dat',status='old')
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

 allocate(m(1:lines),r(1:lines),pres(1:lines))
 do i = 1, rows
  if(trim(header(i))=='mass') m(1:lines) = dat(1:lines,i) * msun
  if(trim(header(i))=='radius') r(1:lines) = dat(1:lines,i) * rsun
  if(trim(header(i))=='pressure') pres(1:lines) = dat(1:lines,i)
 end do

! start setting initial conditions ! -----------------------------------------

!!$ do i = 1, lines
!!$  if(m(i)<PNSmass.and.m(i+1)>=PNSmass)then
!!$   xi1s = (r(i+1)*(PNSmass-m(i))+r(i)*(m(i+1)-PNSmass)) / (m(i+1)-m(i))
!!$   exit
!!$  end if
!!$ end do
!!$ xi1e = r(lines) * 5d1
!!$
!!$!print '(2(1pe13.5e2))',xi1s,xi1e
!!$ call gridset
!!$ call tools
!!$ call metricg
!!$!print *,xi1(ie),xi1(ie+1);stop
!!$!print *,xi1s/rsun,xi1e/rsun;stop

 Ebind = 0d0
 do i = 2, lines
  if(m(i)>1.4d0*msolar)Ebind = Ebind - G*m(i-1)*(m(i)-m(i-1)) / r(i-1)
 end do

 Eexp = Eexp - Ebind

 k = ks
 j = js
 mc(is-1) = PNSmass
 do i = is, ie
  do n = 2, lines-1
   if(xi1(i)>r(n).and.xi1(i)<=r(n+1))then
!    mc(i) = (m(n)*(r(n+1)-xi1(i))+m(n+1)*(xi1(i)-r(n)))/(r(n+1)-r(n))
    mc(i)    = intpol(xi1(i),r(n:n+1),m(n:n+1))
!    mc(i)    = exp(mc(i))
    d(i,j,k) = (mc(i)-mc(i-1)) / dvol(i,j,k)
!    exit
   elseif(xi1(i)>r(lines))then
    d(i,j,k) = 1d-14
    p(i,j,k) = pres(lines) * 1d-3
   end if
   if(x1(i)>(r(n)+r(n-1))*0.5d0.and.x1(i)<=(r(n)+r(n+1))*0.5d0)then
     p(i,j,k) = intpol(log(x1(i)),log((r(n:n+1)+r(n-1:n))*0.5d0),log(pres(n:n+1)))
     p(i,j,k) = exp(p(i,j,k))
   end if
  end do
 end do
!print *,p(1,1,1);stop
 do i = is, is+4
  p(i,j,k) = p(i,j,k) + (gamma-1d0)*Eexp/sum(dvol(is:is+4,js:je,ks))
 end do
 d(is:ie,js:je,k) = spread(d(is:ie,js,k),2,je)
 p(is:ie,js:je,k) = spread(p(is:ie,js,k),2,je)

return

contains

 real*8 function intpol(x,r,u)
  implicit none
  real*8 x
  real*8,dimension(1:2):: r, u

  intpol = ( u(1)*(r(2)-x) + u(2)*(x-r(1)) ) / (r(2)-r(1))
 end function intpol

end subroutine explosion
