!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                    SUBROUTINE ORSZAG-TANG
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set initial condition

subroutine orszagtang

  use grid
  use physval

  implicit none
  
!--------------------------------------------------------------------

! for Orszag-Tang problem : xy surface
  k = ks

  do j = js,je
   do i = is,ie
    d(i,j,k) = gamma*gamma/(4.d0*pi)
    p(i,j,k) = gamma/(4.d0*pi)
    v1(i,j,k)= -dsin(2.d0*pi*x2(j))
    v2(i,j,k)=  dsin(2.d0*pi*x1(i))
    b1(i,j,k)= -dsin(2.d0*pi*x2(j))/dsqrt(4.d0*pi)
    b2(i,j,k)=  dsin(4.d0*pi*x1(i))/dsqrt(4.d0*pi)
   end do
  end do

return
end subroutine orszagtang
