program analytic_shell

! purpose: To calculate the position of shell analytically.

  implicit none

  integer i
  real*8,parameter:: pi = acos(-1d0), year = 36d2*24d0*365.25d0, au = 1.49598073d13
  real*8,parameter:: msun=1.969d33, Mdot0=4d-2*msun/year, vwind0 = 10.d7, tau = 160d0*year
  real*8 theta,deltaM,v0,V,corr, pwind, rshell, Wrot, mshell, mall, rshellold, Mdot, vwind
  logical once

  V = 1d7
  Wrot = 0.9d0
  mall = 0d0
  rshell = 21d3*au
  once = .true.
  open(unit=10,file='homology.dat',status='old')
  read(10,*)
  open(unit=20,file='analytic_homunculus.dat',status='replace')
  write(20,'(4a14)')'theta','rshell','vshell','mshell'
  do
   read(10,*,END=100)theta,deltaM,v0
   Mdot = Mdot0*(1d0-Wrot*sin(theta/180d0*pi)**2d0)
   vwind = vwind0*sqrt(1d0-Wrot*sin(theta/180d0*pi)**2d0)
   pwind = Mdot*vwind*tau! *sqrt(1d0-Wrot*cos(pi/2d0-theta*pi/180d0)**2d0)
   deltaM = deltaM*msun
   corr = 1d99
   do while (abs(corr)>1d-13*abs(V))
!    corr = - (v0*f(V/v0)-pwind*8d0*pi/deltaM) &
!             /fprime(V/v0)
    corr = - (0.5d0*deltaM*v0*f(V/v0)-(1d0-V/vwind)*(pwind-Mdot*tau*V)) &
            /(0.5d0*deltaM*fprime(V/v0)+(Mdot*tau*(1d0-2d0*V/vwind)+pwind/vwind))
    V = V + corr
   end do

   rshellold = rshell
   rshell = V*tau
   mshell = deltaM*0.5d0*fmass(V/v0) + Mdot*tau*(1d0-V/vwind)
   if(theta<90d0)mall = mall + mshell/200d0

   if(once)then
    write(20,'(5(1PE14.6e2))')0d0,rshell/au,V/1e5,mshell/msun,mshell/msun/(1d0+abs(rshellold-rshell)/21d3/au)
    once = .false.
   end if
   write(20,'(5(1PE14.6e2))')theta,rshell/au,V/1e5,mshell/msun,mshell/msun/(1d0+abs(rshellold-rshell)/21d3/au)
  end do

100 close(10)
  close(20)
print *,mall/msun
contains

real*8 function f(x)
implicit none
real*8,intent(in):: x
f = (x*x+4d0*x+6d0)*exp(-x)+2d0*x-6d0
end function f

real*8 function fprime(x)
implicit none
real*8,intent(in)::x
fprime = (-x*x-2d0*x-2d0)*exp(-x)+2d0
end function fprime

real*8 function fmass(x)
implicit none
real*8,intent(in)::x
fmass = -((x*x+2d0*x+2d0)*exp(-x)-2d0)
end function fmass

end program analytic_shell
