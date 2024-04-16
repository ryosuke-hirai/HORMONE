module KHtest_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE KHTESTS
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set initial condition for Kelvini-Helmholtz instability tests.

 subroutine KHtest

  use grid
  use physval

  integer:: i,j,k
  real(8)::rdm
  real(8),parameter:: ptb = 1d-2

!--------------------------------------------------------------------

! for Kelvin-Helmholtz test
  do k = ks-1,ke+1
   do j = js-1,je+1
    do i = is-1,ie+1
     if(abs(x2(j))>0.25d0)then
      d(i,j,k)  = 2d0
      v1(i,j,k) = 0.5d0
     else
      d(i,j,k) = 1d0
      v1(i,j,k) = -0.5d0
     end if

     call random_number(rdm)
     v1(i,j,k) =  v1(i,j,k) + ptb * (2d0*rdm-1d0)
     call random_number(rdm)
     v2(i,j,k) =  ptb * (2d0*rdm-1d0)
     v3(i,j,k) = 0d0
     p(i,j,k)  = 2.5d0
    end do
   end do
  end do

  return
 end subroutine KHtest

end module KHtest_mod
