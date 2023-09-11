module rad_box_mod
 implicit none

contains

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                         SUBROUTINE RAD_BOX
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set initial condition for radiation box.

 subroutine rad_box

  use grid
  use physval

  integer:: i,j,k
  real(8)::rdm
  real(8),parameter:: ptb = 1d-2

!--------------------------------------------------------------------

  do k = ks-1,ke+1
   do j = js-1,je+1
    do i = is-1,ie+1
     d(i,j,k)  = 1d0
     p(i,j,k)  = 1d0
    end do
   end do
  end do

  return
 end subroutine rad_box
 
end module rad_box_mod
