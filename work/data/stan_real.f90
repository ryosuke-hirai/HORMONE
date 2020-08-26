program stan_wind

! purpose: stan_wind

  implicit none

  real*8,parameter:: msun = 1.989d33, year=24d0*36d2*365.25d0, au = 1.49598073d13
  real*8,parameter:: tau=160d0*year, pi = acos(-1d0)
  real*8:: deltaM, Mdot, vwind, u,r,m,v0,t,dt, Wrot, theta, Mdot0, vwind0
  real*8:: ku(1:4),kr(1:4),km(1:4),pall,mall
  logical once

  Mdot0=4d-2*msun/year
  vwind0=1d8
  Wrot = 0.95d0
  once = .true.

  open(unit=10,file='homology.dat',status='old')
  read(10,*)
  open(unit=20,file='stan_homunculus.dat',status='replace')
  write(20,'(4a14)')'theta','rshell','vshell','mshell'

  do
   read(10,*,END=100)theta,deltaM,v0
   Mdot = Mdot0*(1d0-Wrot*sin(theta/180d0*pi)**2d0)
   vwind = vwind0*sqrt(1d0-Wrot*sin(theta/180d0*pi)**2d0)
   deltaM = deltaM*msun
   pall = 0d0;mall = 0d0
   dt = 1d4
   t = 1d3
   m = 1d31!7.340372d-2
   u = 0d0!2.759521d-1
   r = 0d0!2.759521d-3!1d-5

   do while (t<=tau)
    pall = m*u+0.5d0*deltaM*v0*fout(r/t/v0)+Mdot*r
    mall = m + 0.5d0*deltaM*fout2(r/t/v0) +Mdot*r/vwind
!   write(10,'(6(1PE14.6e2))')t,r,u,m,pall,mall
!if(t>1.5d-4)dt=1d-4
    ku(1) = dudt(u,r,m,t)
    kr(1) = drdt(u,r,m,t)
    km(1) = dmdt(u,r,m,t)

    ku(2) = dudt(u+0.5d0*dt*ku(1),r+0.5d0*dt*kr(1),m+0.5d0*dt*km(1),t+0.5d0*dt)
    kr(2) = drdt(u+0.5d0*dt*ku(1),r+0.5d0*dt*kr(1),m+0.5d0*dt*km(1),t+0.5d0*dt)
    km(2) = dmdt(u+0.5d0*dt*ku(1),r+0.5d0*dt*kr(1),m+0.5d0*dt*km(1),t+0.5d0*dt)

    ku(3) = dudt(u+0.5d0*dt*ku(2),r+0.5d0*dt*kr(2),m+0.5d0*dt*km(2),t+0.5d0*dt)
    kr(3) = drdt(u+0.5d0*dt*ku(2),r+0.5d0*dt*kr(2),m+0.5d0*dt*km(2),t+0.5d0*dt)
    km(3) = dmdt(u+0.5d0*dt*ku(2),r+0.5d0*dt*kr(2),m+0.5d0*dt*km(2),t+0.5d0*dt)

    ku(4) = dudt(u+dt*ku(3),r+dt*kr(3),m+dt*km(3),t+dt)
    kr(4) = drdt(u+dt*ku(3),r+dt*kr(3),m+dt*km(3),t+dt)
    km(4) = dmdt(u+dt*ku(3),r+dt*kr(3),m+dt*km(3),t+dt)
!print *,ku,kr,km
!stop
    u = u + dt*(ku(1)+2d0*ku(2)+2d0*ku(3)+ku(4))/6d0
    r = r + dt*(kr(1)+2d0*kr(2)+2d0*kr(3)+kr(4))/6d0
    m = m + dt*(km(1)+2d0*km(2)+2d0*km(3)+km(4))/6d0

    t = t + dt

   end do
   if(once)then
    write(20,'(5(1PE14.6e2))')0d0,r/au,u/1e5,m/msun
    once = .false.
   end if
   write(20,'(5(1PE14.6e2))')theta,r/au,u/1e5,m/msun
  end do

100 close(10)
   close(20)

contains

real*8 function dmdt(uu,rr,mm,tt)
 implicit none

 real*8,intent(in):: uu,rr,mm,tt

 dmdt = Mdot/vwind*(vwind-uu)+deltaM*rr*rr/(2d0*(v0*tt)**3d0)*exp(-rr/(v0*tt))*(uu-rr/tt)

end function dmdt

real*8 function dudt(uu,rr,mm,tt)
 implicit none

 real*8,intent(in):: uu,rr,mm,tt

 dudt = (Mdot/vwind*(vwind-uu)**2d0-deltaM*rr*rr/(2d0*(v0*tt)**3d0)*exp(-rr/(v0*tt))*(uu-rr/tt)**2d0)/mm! + uu*(Mdot/vwind*(vwind-uu)+deltaM*rr*rr/(2d0*(v0*tt)**3d0)*exp(-rr/(v0*tt))*(uu-rr/tt))

end function dudt

real*8 function drdt(uu,rr,mm,tt)
 implicit none

 real*8,intent(in):: uu,rr,mm,tt

 drdt = uu

end function drdt

real*8 function fout(x)
implicit none
real*8,intent(in)::x
fout = (x*x*x+3d0*x*x+6d0*x+6d0)*exp(-x)
end function fout

real*8 function fout2(x)
implicit none
real*8,intent(in)::x
fout2 = (x*x+2d0*x+2d0)*exp(-x)
end function fout2

end program stan_wind
