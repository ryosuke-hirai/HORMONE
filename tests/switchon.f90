!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE SWITCHON
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set initial condition

subroutine switchon

  use grid
  use physval

  implicit none
  
!--------------------------------------------------------------------

! for switch-on shock
k = ks ; j = js
 do i = is,in/2
  d(i,j,k) = 1.d0
  v1(i,j,k) = 0.d0
  v2(i,j,k) = 0.d0
  v3(i,j,k) = 0.d0
  b1(i,j,k) = 3.d0
  b2(i,j,k) = 1.d0
  b3(i,j,k) = 0.d0
  p(i,j,k)  = 1.d0
 end do
 do i = in/2+1,ie
  d(i,j,k) = 2.62282566928d0
  v1(i,j,k) = -2.19684296835d0
  v2(i,j,k) = -1.57158410509d0
  v3(i,j,k) = 0.d0
  b1(i,j,k) = 3.d0
  b2(i,j,k) = -0.86d0
  b3(i,j,k) = 0.d0
  p(i,j,k)  = 8.93021765330d0
 end do

return
end subroutine switchon
