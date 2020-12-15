!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!
!                           SUBROUTINE CONSERVE
!
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

! PURPOSE: To convert physical values to conserved values

subroutine conserve

  use grid
  use physval

  implicit none

!---------------------------------------------------------------------------

  ! set U
!$omp parallel do private(i,j,k)
  do k = ks,ke
   do j = js,je
    do i = is,ie
     u(i,j,k,1) = d(i,j,k)
     u(i,j,k,2) = d(i,j,k) * v1(i,j,k)
     u(i,j,k,3) = d(i,j,k) * v2(i,j,k)
     u(i,j,k,4) = d(i,j,k) * v3(i,j,k)
     u(i,j,k,5) = b1(i,j,k)
     u(i,j,k,6) = b2(i,j,k)
     u(i,j,k,7) = b3(i,j,k)
     u(i,j,k,8) = e(i,j,k)
     u(i,j,k,9) = phi(i,j,k)
    end do
   end do
  end do
!$omp end parallel do

return
end subroutine conserve


