!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE INITIALCONDITION
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set initial condition

subroutine initialcondition

 use funcs
 use settings,only:start
 use grid
 use physval
 use pressure_mod

 use constants
 use gravmod
 use merger_mod

 implicit none

 real*8:: Mdot, vinf

!----------------------------------------------------------------------------

 if(start==0)then

!  call KHtest

!  call shocktube

!  call orszagtang
  
!  call switchon

!  call sedovtaylor

!  call explosion

!  call companion

!  call headoncollision

!  call spinupRSG

  call failedSN

! to set internal energy consistently -------------------------------------- !
  call meanmolweight                                                         !
  do k=ks,ke ; do j=js,je ; do i=is,ie                                       !
   e(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k)) &                 !
  + d(i,j,k)*0.5d0*( pw(2,v1(i,j,k)) + pw(2,v2(i,j,k)) + pw(2,v3(i,j,k)) ) & !
  +          0.5d0*( pw(2,b1(i,j,k)) + pw(2,b2(i,j,k)) + pw(2,b3(i,j,k)) )   !
  end do ; end do ; end do                                                   !
! -------------------------------------------------------------------------- !

 else

!  call failedSN
  call restart
  call meanmolweight
xi1e = xi1(ie)
! to get rid of outside shell
  if(time==inifile)then
   Mdot = 1d-6*msun/year
   vinf = 1.5d7
   do k = ks, ke
    do j = js, je
     do i = is, ie
      if(d(i,j,k)*v3(i,j,k)<3d-5)then
       d (i,j,k) = Mdot/(4d0*pi*x1(i)*x1(i)*vinf)
       v1(i,j,k) = vinf
       v2(i,j,k) = 0d0
       v3(i,j,k) = 0d0
       p (i,j,k) = -d(i,j,k)*grvphi(i,j,k)
       spc(1:8,i,j,k) = spc(1:8,ie,je,ks)
       imu(i,j,k) = 0.25d0*(6d0*spc(1,i,j,k)+spc(2,i,j,k)+2d0)
       eint(i,j,k) = eos_e(d(i,j,k),p(i,j,k),T(i,j,k),imu(i,j,k))
       e (i,j,k) = eint(i,j,k) + 0.5d0*d(i,j,k)*v1(i,j,k)*v1(i,j,k)
      else
       eint(i,j,k) = e(i,j,k) - 0.5d0*d(i,j,k)*(v1(i,j,k)**2d0+v2(i,j,k)**2d0+v3(i,j,k)**2d0)
       p (i,j,k) = eos_p(d(i,j,k),eint(i,j,k),T(i,j,k),imu(i,j,k))
       e (i,j,k) = e(i,j,k) - 0.5d0*d(i,j,k)*(v1(i,j,k)**2d0+v2(i,j,k)**2d0)
       v1(i,j,k) = 0d0
       v2(i,j,k) = 0d0
      end if
     end do
    end do
   end do
  end if

 end if

return
end subroutine initialcondition
