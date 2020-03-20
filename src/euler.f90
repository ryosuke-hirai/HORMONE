!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                        SUBROUTINE EULER
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To integrate numflux using Euler method

subroutine euler

  use grid
  use physval

  implicit none

!--------------------------------------------------------------------

  do ufn = 1,9
   do k = ks,ke
    do j = js,je
     do i = is,ie
      u(i,j,k,ufn) = u(i,j,k,ufn) - dt * &
           ( idetg1(i) * &
             (detg1(i  )*flux1(i,j,k,ufn)-detg1(i-1  )*flux1(i-1,j,k,ufn)) &
           + idetg2(i,j) * &
             (detg2(i,j)*flux2(i,j,k,ufn)-detg2(i,j-1)*flux2(i,j-1,k,ufn)) &
           + idetg3(i,j,k) * (flux3(i,j,k,ufn)-flux3(i,j,k-1,ufn)) &
           + src(i,j,k,ufn) )
     end do
    end do
   end do
  end do

!  call primitive

return
end subroutine euler
