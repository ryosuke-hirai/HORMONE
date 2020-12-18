program lightcurve

! purpose: To calculate a simple light curve from the hydro data.

 use grid
 use physval
 use constants
 use readbin_mod
 
 implicit none

 character*50::binfile,outfile
 real*8:: rad, lum
 integer:: angle

!-----------------------------------------------------------------------------

 angle = 60
 rad = dble(angle)*pi/180d0

 call allocate_readgrid
 
 write(outfile,'(i2.2,"deg_lightcurve.dat")')angle
 open(unit=15,file=outfile,status='replace')

 do n = 300, 1600, 10
  write(binfile,'("bin",i11.11,"hr.dat")')n
  
  call readbin(binfile)
  call get_luminosity(rad,lum)
  print*,n,lum
  write(15,'(i5,1PE14.6e2)')n,lum
  flush(15)
 end do

 close(15)
 
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
end module utils

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                      SUBROUTINE GET_LUMINOSITY
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
 real*8:: tau, tau0, dr, dtheta, dA, rho, XX, TT, dzabs, lumt
 real*8,dimension(1:3):: dz=0d0, x, xp
 integer:: ii, jj, kk, iin, jjn, kkn

!-----------------------------------------------------------------------------

 tau0 = 1d0 ! optical depth at photosphere
 iin = 200 ! radial resolution
 jjn = 100 ! angular resolution
 kkn = 1000 ! depth resolution
 dr = xi1(ie)/dble(iin)
 dtheta = 2d0*pi/dble(jjn)
 dz = (/0d0,0d0,-xi1(ie)/dble(kkn)/)
 dz = rotx(dz,angle)
 dzabs = norm2(dz)
 lum = 0d0
 lumt = 0d0

!$omp parallel do default(none) reduction(+:lumt) &
!$omp private(ii,jj,kk,dA,x,xp,tau,rho,XX,TT) &
!$omp shared(dr,dtheta,dzabs,dz,tau0,angle,iin,jjn,kkn,xi1,ie)
 do jj = 1, jjn
  do ii = 1, iin
! set initial position
   tau = 0d0
   x(1) = dr*(dble(ii)-0.5d0)*cos(jj*dtheta)
   x(2) = dr*(dble(ii)-0.5d0)*sin(jj*dtheta)
   x(3) = xi1(ie)
   dA = (dble(ii)-0.5d0)*dr**2*dtheta
   x = rotx(x,angle)
   xp = carpol(x)
! find photosphere
   do kk = 1, kkn
    call get_local_val(xp,rho,XX,TT)
    tau = tau + rho*kap_es(XX)*dzabs
 !   print*,ii,jj,tau,xp(1)/au
    if(tau>tau0)exit
    x = x + dz
    xp = carpol(x)
   end do
! add black body flux
   lumt = lumt+dA*sigma*TT**4
  end do
 end do
!$omp end parallel do

 lum = lumt
 
 return
end subroutine get_luminosity

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!                      SUBROUTINE GET_LOCAL_VAL
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: Get local values of density, temperature and hydrogen fraction

 subroutine get_local_val(xp,rho,XX,TT)

  use grid
  use physval
  use constants,only:pi

  implicit none

  real*8,intent(in):: xp(1:3)
  real*8,intent(out):: rho, XX, TT
  real*8:: xp2
  integer:: ii,jj,kk

 !-----------------------------------------------------------------------------

  if(xp(1)>=xi1(ie))then
   rho = 0d0
   XX  = 0.7d0
   TT  = 1d3
   return
  end if
  kk = ks
  xp2 = min(xp(2),pi-xp(2))
  do jj = js, je
   if(xi2(jj)>=xp2)exit
  end do
  do ii = is, ie
   if(xi1(ii)>=xp(1))exit
  end do
  rho = d(ii,jj,kk)
  XX  = spc(1,ii,jj,kk)
  TT  = T(ii,jj,kk)
  
 return
 end subroutine get_local_val


