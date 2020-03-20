!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                            SUBROUTINE COMPANION
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To input the initial condition for 2D companion star.

subroutine companion

 use grid
 use physval
 use constants
 use gravmod,only:mc
 use ejmod

 implicit none

 integer lines, rows, n, nn
 character*10000 dumc
 character*20,allocatable:: header(:),dum(:)
 real*8,allocatable,dimension(:,:):: dat
 real*8,allocatable,dimension(:):: m, r, pres
 real*8,allocatable,dimension(:,:,:):: pres0, rho0, b10, b20
 real*8 rstar

!-----------------------------------------------------------------------------

! reading data from datafile ! -----------------------------------------------
 open(unit=40,file='16-10_025_star2.dat',status='old')
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

 rstar = r(lines)
 k = ks
 j = js
 mc(is-1) = 0d0
 do i = is, ie
  do n = 2, lines-1
   if(xi1(i)>r(n).and.xi1(i)<=r(n+1))then
!    mc(i) = (m(n)*(r(n+1)-xi1(i))+m(n+1)*(xi1(i)-r(n)))/(r(n+1)-r(n))
    mc(i)    = intpol(xi1(i),r(n:n+1),m(n:n+1))
!    mc(i)    = exp(mc(i))
    d(i,j,k) = (mc(i)-mc(i-1)) / sum(dvol(i,js:je,k))
!    exit
   elseif(xi1(i)>r(lines))then
    d(i,j,k) = 1d-17
    p(i,j,k) = pres(lines) * 1d-04
   end if
   if(x1(i)>(r(n)+r(n-1))*0.5d0.and.x1(i)<=(r(n)+r(n+1))*0.5d0)then
     p(i,j,k) = intpol(log(x1(i)),log((r(n:n+1)+r(n-1:n))*0.5d0),log(pres(n:n+1)))
     p(i,j,k) = exp(p(i,j,k))
   end if
  end do
 end do

 d(is:ie,js:je,k) = spread(d(is:ie,js,k),2,je)
 p(is:ie,js:je,k) = spread(p(is:ie,js,k),2,je)

!!$ call readejecta


!!$allocate( pres0(ie,je,ks), rho0(ie,je,ks), b10(ie,je,ks), b20(ie,je,ks) )
!!$k = ks
!!$  do i = is,ie
!!$   do j = js,je
!!$    pres0(i,j,k) = sqrt( sep**2.+x1(i)**2. -2.d0*x1(i)*sep*cosc(j) )
!!$    rho0(i,j,k)  = ejectadistance / pres0(i,j,k)
!!$    b10(i,j,k)   = sinc(j) * sep / pres0(i,j,k)
!!$    b20(i,j,k)   = sqrt(1.d0-b10(i,j,k)**2.d0)
!!$   end do
!!$  end do
!!$
!!$  time = tstart / ejectadistance * (sep-rstar*1.01)
!!$
!!$  do i = is,ie
!!$   do j = js,je/2
!!$
!!$    if(time*rho0(i,j,k)>=tstart)then
!!$     do nn = 1, count-1
!!$      if(time*rho0(i,j,k)>=t_ej(nn).and.time*rho0(i,j,k)<t_ej(nn+1))then
!!$       d(i,j,k)   = (d_ej(nn)   * (t_ej(nn+1)-time*rho0(i,j,k))  &
!!$                   - d_ej(nn+1) * (t_ej(nn)  -time*rho0(i,j,k))) &
!!$                   / (t_ej(nn+1)-t_ej(nn)) &
!!$                   * rho0(i,j,k)**3.d0
!!$       p(i,j,k)   = (p_ej(nn)   * (t_ej(nn+1)-time*rho0(i,j,k))  &
!!$                   - p_ej(nn+1) * (t_ej(nn)  -time*rho0(i,j,k))) &
!!$                   / (t_ej(nn+1)-t_ej(nn)) &
!!$                   * rho0(i,j,k)**(3.d0*gamma)
!!$       e(i,j,k)   = (e_ej(nn)   * (t_ej(nn+1)-time*rho0(i,j,k))  &
!!$                   - e_ej(nn+1) * (t_ej(nn)  -time*rho0(i,j,k))) &
!!$                   / (t_ej(nn+1)-t_ej(nn)) &
!!$                   * rho0(i,j,k)**5.d0
!!$       v1(i,j,k) =  (v_ej(nn)   * (t_ej(nn+1)-time*rho0(i,j,k))  &
!!$                   - v_ej(nn+1) * (t_ej(nn)  -time*rho0(i,j,k))) &
!!$                   / (t_ej(nn+1)-t_ej(nn)) &
!!$                   * (-b20(i,j,k))! * rho0(i,j,k)
!!$       v2(i,j,k) =  (v_ej(nn)   * (t_ej(nn+1)-time*rho0(i,j,k))  &
!!$                   - v_ej(nn+1) * (t_ej(nn)  -time*rho0(i,j,k))) &
!!$                   / (t_ej(nn+1)-t_ej(nn)) &
!!$                   * b10(i,j,k)! * rho0(i,j,k)
!!$       if(p_ej(nn)>pmax*1.d-3)then
!!$        p(i,j,k) = pmax !for shock condition
!!$       end if
!!$      end if
!!$
!!$     end do
!!$    end if
!!$
!!$   end do
!!$  end do

return

contains

 real*8 function intpol(x,r,u)
  implicit none
  real*8 x
  real*8,dimension(1:2):: r, u

  intpol = ( u(1)*(r(2)-x) + u(2)*(x-r(1)) ) / (r(2)-r(1))
 end function intpol

end subroutine companion
