!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                    SUBROUTINE ORSZAG-TANG
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To set initial condition

subroutine orszagtang

  use grid
  use physval
  use constants,only:pi

  implicit none
  
!--------------------------------------------------------------------

! for Orszag-Tang problem : xy surface
  k = ks

  do j = js,je
   do i = is,ie
    d(i,j,k) = gamma**2/(4.d0*pi)
    p(i,j,k) = gamma/(4.d0*pi)
    v1(i,j,k)= -sin(2.d0*pi*x2(j))
    v2(i,j,k)=  sin(2.d0*pi*x1(i))
    v3(i,j,k)= 0d0
    b1(i,j,k)= -sin(2.d0*pi*x2(j))/sqrt(4.d0*pi)
    b2(i,j,k)=  sin(4.d0*pi*x1(i))/sqrt(4.d0*pi)
    b3(i,j,k)= 0d0
   end do
  end do

return
end subroutine orszagtang
