module KHinstability_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE KHINSTABILITY
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set initial condition for Kelvin-Helmholtz instability tests.

 subroutine KHinstability

  use constants,only:pi
  use grid
  use physval

  integer:: i,j,k
  real(8):: rho1,rho2,rhom,vv1,vv2,vvm
  real(8),parameter:: ptb = 1d-2, width=0.025d0

!--------------------------------------------------------------------

  species(1) = 'scalar1'
  species(2) = 'scalar2'
  rho1 = 1d0
  rho2 = 2d0
  rhom = 0.5d0*(rho1-rho2)
  vv1  = 0.5d0
  vv2  =-0.5d0
  vvm  = 0.5d0*(vv1-vv2)

! for Kelvin-Helmholtz test
  do k = ks-1,ke+1
   do j = js-1,je+1
    do i = is-1,ie+1
     if(abs(x2(j))>0.25d0)then
      v1(i,j,k) = vv1 - vvm *exp((-abs(x2(j))+0.25d0)/width)
      d(i,j,k) = rho1 - rhom*exp((-abs(x2(j))+0.25d0)/width)
      spc(1,i,j,k) = 1d0
      spc(2,i,j,k) = 0d0
     else
      v1(i,j,k) = vv2 + vvm *exp((abs(x2(j))-0.25d0)/width)
      d(i,j,k) = rho2 + rhom*exp((abs(x2(j))-0.25d0)/width)
      spc(1,i,j,k) = 0d0
      spc(2,i,j,k) = 1d0
     end if

     v2(i,j,k) = ptb * sin(4d0*pi*x1(i))
     v3(i,j,k) = 0d0
     p (i,j,k) = 2.5d0
    end do
   end do
  end do

  return
 end subroutine KHinstability

end module KHinstability_mod
