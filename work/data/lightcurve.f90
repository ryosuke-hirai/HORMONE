program lightcurve

! purpose: To calculate a simple light curve from the hydro data.

 use grid
 use physval
 use constants
 use readbin_mod
 
 implicit none

 character*50::binfile,outfile
 real*8:: rad, lum
 integer,parameter:: iangle=0, fangle=0, dangle=10
 integer:: angle, unitn(0:fangle)

!-----------------------------------------------------------------------------

 call allocate_readgrid
 
 do angle = iangle, fangle, dangle

  write(outfile,'(i2.2,"deg_lightcurve.dat")')angle
  open(newunit=unitn(angle),file=outfile,status='replace')

 end do

 do n = 0, 500, 1
  write(binfile,'("bin",i11.11,"hr.dat")')n
  call readbin(binfile)
  
  do angle = iangle, fangle, dangle
   rad = dble(angle)*pi/180d0
   call get_luminosity(rad,lum)
   print*,n,angle,lum
   write(unitn(angle),'(i5,1PE14.6e2)')n,lum
   flush(unitn(angle))
  end do
  
 end do

 do angle = iangle, fangle, dangle
  close(unitn(angle))
 end do

end program lightcurve


module utils
  implicit none

 contains

! electron scattering opacity
 pure real*8 function kap_es(x)
  implicit none
  real*8,intent(in):: x
  kap_es = 0.2d0*(1d0+x)
 end function kap_es

 ! rotate about the x axis
 pure function rotx(x,theta) result(xp)
  implicit none
  real*8,intent(in):: x(1:3), theta
  real*8:: xp(1:3)
  xp(1) = x(1)
  xp(2) = cos(theta)*x(2) - sin(theta)*x(3)
  xp(3) = sin(theta)*x(2) + cos(theta)*x(3)
 end function rotx

! rotate about the y axis
 pure function roty(x,theta) result(xp)
  implicit none
  real*8,intent(in):: x(1:3), theta
  real*8:: xp(1:3)
  xp(1) = cos(theta)*x(1) + sin(theta)*x(3)
  xp(2) = x(2)
  xp(3) =-sin(theta)*x(1) + cos(theta)*x(3)
 end function roty

! rotate about the z axis
 pure function rotz(x,theta) result(xp)
  implicit none
  real*8,intent(in):: x(1:3), theta
  real*8:: xp(1:3)
  xp(1) = cos(theta)*x(1) - sin(theta)*x(2)
  xp(2) = sin(theta)*x(1) + cos(theta)*x(2)
  xp(3) = x(3)
 end function rotz

! convert polar to cartesian coordinates
 pure function polcar(xp) result(x)
  implicit none
  real*8,intent(in):: xp(1:3)
  real*8:: x(1:3)
  x(1) = xp(1)*sin(xp(2))*cos(xp(3))
  x(2) = xp(1)*sin(xp(2))*sin(xp(3))
  x(3) = xp(1)*cos(xp(2))
 end function polcar

! convert cartesian to polar coordinates
 pure function carpol(x) result(xp)
  implicit none
  real*8,intent(in):: x(1:3)
  real*8:: xp(1:3)
  xp(1) = norm2(x)
  xp(2) = acos(x(3)/xp(1))
  xp(3) = atan2(x(2),x(1))
 end function carpol

! linear interpolation
 pure function intpol(x,y,z) result(val)
  implicit none
  real*8,intent(in):: x(1:2), y(1:2)
  real*8,intent(in):: z
  real*8:: val
  val = ((x(2)-z)*y(1)+(z-x(1))*y(2))/(x(2)-x(1))
 end function intpol

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                      SUBROUTINE GEOMETRICAL_SERIES
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate dxi's in a geometrical series.

subroutine geometrical_series(dxi,xmin,is,ie,xis,xie)

 implicit none

 integer,intent(in):: is,ie
 real*8,intent(in):: xis,xie,xmin
 real*8,intent(inout),allocatable:: dxi(:)
 integer i
 real*8 xrng, irng, xr, xrnew, xrmax, err, maxerr, fx, dfx

!-----------------------------------------------------------------------------

 xrmax = 1.15d0
 maxerr = 1d-10

 xr = 1.01d0
 xrng = xie - xis ; irng = dble(ie - is + 1)

 if(xrng/irng<xmin)then
  print *,"Error from geometrical_series ;"
  print *,"xmin should be smaller or uniform mesh should be chosen",xmin
  stop
 end if

 do i = 1, 10000000
  fx = (xr-1d0)*xrng - xmin * (xr**irng-1d0)
  dfx = xrng - irng * xmin * xr**(irng-1d0)

  xrnew = xr - fx/dfx

  if(abs((xrnew-xr)/xr)<maxerr)then
   xr = xrnew ; exit
  end if
  if(xrnew<1d0)xrnew = 2d0

  xr = xrnew
 end do

 if(xr>xrmax)then
  print *,"xmin too small", xmin, xr
  stop
 end if

 dxi(is) = xmin
 do i = is+1, ie
  dxi(i) = dxi(i-1) * xr
 end do
 dxi(is-1) = dxi(is) ; dxi(is-2) = dxi(is+1)
 dxi(ie+1) = dxi(ie)*xr ; dxi(ie+2) = dxi(ie)*xr*xr

 if(xr-1d0<maxerr) dxi = (xie-xis) / irng

return
end subroutine geometrical_series

end module utils

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                       SUBROUTINE GET_LUMINOSITY
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To calculate luminosity assuming black body radiation

subroutine get_luminosity(angle,lum)

 use grid
 use physval
 use constants
 use utils

 implicit none

 real*8,intent(in):: angle
 real*8,intent(out):: lum
 real*8:: tau, tau0, drmin, dtheta, dA, rho, XX, TT, dzabs, lumt, dzabs0
 real*8,dimension(1:3):: dz, x, xp, dz0
 real*8,parameter:: tauerr=1d-3
 real*8,allocatable:: dr(:), rr(:)
 integer:: ii, jj, kk, iin, jjn, kkn

!-----------------------------------------------------------------------------

 tau0 = 1d0 ! optical depth at photosphere
 iin = 300 ! radial resolution
 jjn = 10 ! angular resolution
 kkn = 100 ! depth resolution
 allocate(dr(-1:iin+2),rr(-1:iin+2))
 drmin = dxi1(is)*0.2d0
 call geometrical_series(dr,drmin,1,iin,0d0,xi1(ie))
 rr(0) = 0d0
 do ii = 1, iin
  rr(ii) = rr(ii-1)+dr(ii)
 end do
 dtheta = 2d0*pi/dble(jjn)
 dz0 = (/0d0,0d0,-xi1(ie)/dble(kkn)/)
 dz0 = rotx(dz0,angle)
 dzabs0 = norm2(dz0)
 lum = 0d0
 lumt = 0d0

!$omp parallel do default(none) reduction(+:lumt) &
!$omp private(ii,jj,kk,dA,x,xp,tau,rho,XX,TT,dz,dzabs) &
!$omp shared(dtheta,tau0,angle,iin,jjn,kkn,xi1,ie,dz0,dzabs0,rr)
 do jj = 1, jjn
  do ii = 1, iin
! set initial position
   tau = 0d0;dz=dz0;dzabs=dzabs0
   x(1) = 0.5d0*(rr(ii-1)+rr(ii))*cos(jj*dtheta)
   x(2) = 0.5d0*(rr(ii-1)+rr(ii))*sin(jj*dtheta)
   x(3) = sqrt(xi1(ie)**2-x(1)**2-x(2)**2)
   dA = 0.5d0*(rr(ii)**2-rr(ii-1)**2)*dtheta
   x = rotx(x,angle)
   xp = carpol(x)
! find photosphere
   do while (dot_product(x,dz0)<xi1(ie)*dzabs0)
    x = x + dz
    xp = carpol(x)
    call get_local_val(xp,rho,XX,TT)
    tau = tau + rho*kap_es(XX)*dzabs

    if(tau-tau0>tauerr*tau0)then ! tan ran over tau0
     tau = tau - rho*kap_es(XX)*dzabs
     x = x - dz
     xp = carpol(x)
     dz = dz*0.1d0
     dzabs = norm2(dz)!dzabs*0.1d0
     cycle
    elseif(tau>tau0.and.abs(tau-tau0)<tauerr*tau0)then ! tau=tau0
! add black body flux
     call get_local_val(xp,rho,XX,TT)
     lumt = lumt+dA*sigma*TT**4
!     if(TT>0d0)print'(2i5,3(1PE13.5e2))',ii,jj,x(3),TT,tau
     exit
    end if
   end do

  end do
 end do
!$omp end parallel do

 lum = lumt

 return
end subroutine get_luminosity

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                        SUBROUTINE GET_LOCAL_VAL
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Get local values of density, temperature and hydrogen fraction

 subroutine get_local_val(xp,rho,XX,TT)

  use grid
  use physval
  use constants,only:pi
  use utils,only:intpol
  
  implicit none

  real*8,intent(in):: xp(1:3)
  real*8,intent(out):: rho, XX, TT
  real*8:: xp1, xp2, rhot(1:2), XXt(1:2), TTt(1:2)
  integer:: ii,jj,kk

 !-----------------------------------------------------------------------------

  if(xp(1)>=x1(ie))then
   rho = 0d0
   XX  = 0.7d0
   TT  = 1d-30
   return
  end if
  kk = ks
  xp1 = xp(1)
  xp2 = min(xp(2),pi-xp(2))
  do jj = js, je
   if(x2(jj)>=xp2)exit
  end do
  if(jj==js)xp2=x2(js)
  do ii = is, ie
   if(x1(ii)>=xp1)then
    if(ii==is)xp1 = x1(ii)
    rhot(1) = intpol(x1(ii-1:ii),d(ii-1:ii,jj-1,kk),xp1)
    rhot(2) = intpol(x1(ii-1:ii),d(ii-1:ii,jj  ,kk),xp1)
    rho = intpol(x2(jj-1:jj),rhot(1:2),xp2)
    XXt(1) = intpol(x1(ii-1:ii),spc(1,ii-1:ii,jj-1,kk),xp1)
    XXt(2) = intpol(x1(ii-1:ii),spc(1,ii-1:ii,jj  ,kk),xp1)
    XX = intpol(x2(jj-1:jj),XXt(1:2),xp2)
    TTt(1) = intpol(x1(ii-1:ii),T(ii-1:ii,jj-1,kk),xp1)
    TTt(2) = intpol(x1(ii-1:ii),T(ii-1:ii,jj  ,kk),xp1)
    TT = intpol(x2(jj-1:jj),TTt(1:2),xp2)
    exit
   end if
  end do

 return
 end subroutine get_local_val


